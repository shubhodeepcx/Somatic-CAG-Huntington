"""
hdmodel.py -- Core multiscale model of somatic CAG-repeat expansion in
Huntington's disease striatal projection neurons (SPNs).

Layers
------
L1  Biochemical : length-dependent per-cell expansion rate lambda(L) built from
                  (i) MutSbeta site multiplicity n(L) and (ii) the equilibrium
                  occupancy of extrahelical 1-unit slip-outs phi(L).
L2  Cellular    : exact stochastic (Gillespie) simulation of L(t) in single SPNs,
                  with a length-dependent elimination hazard.
L3  PK/PD       : MSH3 protein lowering -> fractional suppression of lambda,
                  using the published mouse di-siRNA dose-response regression.

All literature-derived constants are collected in `LIT` with their source tags so
that every number used in the model is traceable. Parameters that are *fitted*
in this work are returned by `calibrate()`; parameters that are *assumed* are
flagged `ASSUMED` and are all covered by the global sensitivity analysis.

No experimental data are redistributed here: calibration uses published summary
statistics only (see paper, Table 1).
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Literature-derived constants
# ---------------------------------------------------------------------------

LIT = {
    # --- Handsaker et al., Cell 2025 (PMID 39824182) : SPN expansion kinetics
    "t_40_60_years": 50.0,      # mean time for an SPN to expand 40 -> 60 CAG
    "t_60_80_years": 12.0,      # mean time to expand 60 -> 80 CAG
    "t_80_150_years": 4.0,      # "just a few years"; taken as 4 y, swept 2-8 y
    "rate_rise_70_90_min": 6.0, # mutation rate rises >= 6-fold across 70-90 CAG
    "L_phaseC": 150.0,          # onset of cell-autonomous transcriptional change
    "frac_life_below_150": 0.964,  # 96.4% (SD 2.0) of SPN life spent at 42-150
    "derepress_frac_150_250": 0.04,   # de-repression rare (<4%) at 150-250
    "derepress_frac_350plus": 0.50,   # very common (>50%) at >=350
    "max_observed_CAG": 842.0,

    # --- Li, Wang & Yang, bioRxiv 2023 (PMID 38168405; PREPRINT) : MutSbeta
    "mutsb_footprint_repeats": 17.0,   # ~17 CNG repeats occluded per MutSbeta at a stem
    "mutsb_Kd_stem_pM": 129.0,         # (CAG)18 hairpin stem
    "mutsb_Kd_end_pM": 190.0,          # (CAG)8 hairpin end
    "mutsb_1us_fraction": 0.90,        # >90% of stem-bound particles engage a 1-unit slip-out

    # --- Senoussi et al., Nat Commun 2025;16:9445 : expansion step size
    "step_repeats": 4.0,               # (CTG)4 = 12 nt insertion observed in vitro

    # --- O'Reilly et al. / Sci Rep 2024, doi:10.1038/s41598-024-52667-3 : MSH3 di-siRNA dose-response
    #     linear regression of % prevention of striatal expansion against
    #     striatal MSH3 protein (% of vehicle): y-intercept 92.4, x-intercept 110.4
    "pd_y_intercept": 92.4,
    "pd_x_intercept": 110.4,
    "pd_r2": 0.65,
    "msh3_protein_10nmol_pct": 25.0,   # striatal MSH3 protein, % of vehicle, 10 nmol ICV
    "msh3_protein_3nmol_pct": 59.0,    # 3 nmol ICV
    "expansion_prevented_10nmol": 78.1,
    "expansion_prevented_3nmol": 45.1,

    # --- Frei et al., Nucleic Acid Ther 2025 : intrathecal ASO in NHP striatum
    "nhp_putamen_spn_residual": 0.71,  # residual target RNA in putamen SPNs
    "nhp_caudate_residual_snrna": 0.57,
    "nhp_caudate_residual_qpcr": 0.45,
    "nhp_putamen_residual_bulk": 0.75,

    # --- Johnsen et al., Theranostics 2018;8:3416 : TfR-targeted AuNP, 8 h post-IV
    "aunp_parenchyma_pctID_g": {"bivalent_highaff": 0.04,
                                "bivalent_lowaff": 0.08,
                                "monovalent": 0.23},
    "aunp_capillary_pctID_g": {"bivalent_highaff": 0.2,
                               "bivalent_lowaff": 0.5,
                               "monovalent": 1.1},
    "aunp_hydrodynamic_nm": 75.0,

    # --- Orth & Schwenke, PLOS Curr Huntingt Dis 2011 (Ranen analogous model)
    #     AAO = 90.3918 + 0.3293*parental AAO - 1.3996*CAG
    "ranen_intercept": 90.3918,
    "ranen_parent_coef": 0.3293,
    "ranen_cag_coef": -1.3996,
    "ranen_mean_parent_aao": 45.0,   # ASSUMED cohort mean; swept in 04_reconcile.py

    # --- Enroll-HD / Warner et al. : CAP score
    "cap_offset": 33.66,
}

# ---------------------------------------------------------------------------
# Assumed parameters (all swept in the global sensitivity analysis)
# ---------------------------------------------------------------------------

ASSUMED = {
    "L50_ceiling": 400.0,  # length at which slip-out occupancy saturates (soft ceiling)
    "w_switch": 6.0,       # width (repeats) of the isolated->clustered slip-out switch
    "L_death_50": 300.0,   # CAG length at which elimination hazard is half-maximal
    "death_width": 40.0,   # steepness (CAG) of the elimination hazard sigmoid
    "h_death_max": 6.0,    # maximal elimination hazard (per year)
    "contraction_ratio": 0.10,  # contraction events per expansion event
    "contraction_step": 2.0,    # mean repeats lost per contraction event
    "L_max": 1000.0,       # numerical ceiling on tract length
}


# ---------------------------------------------------------------------------
# Layer 1 : biochemical rate law
# ---------------------------------------------------------------------------

def n_sites(L, footprint=None):
    """Number of MutSbeta molecules a tract of length L can accommodate in tandem.

    Li et al. report a stem footprint of ~17 CNG repeats per MutSbeta, and
    EMSA stoichiometries of 1-2 molecules on (CAG)25 and up to 4 on (CAG)60,
    i.e. roughly L/17 rounded down. We use the continuous form L/footprint so
    that the rate law is smooth, with a floor of 1.
    """
    w = LIT["mutsb_footprint_repeats"] if footprint is None else footprint
    return np.maximum(L / w, 1.0)


def phi_slipout(L, beta, L50):
    """Equilibrium probability that a register is in the extrahelical 1-unit
    slip-out state that MutSbeta recognises.

    Modelled as a two-state equilibrium whose free-energy difference falls
    linearly with tract length (each additional repeat adds stacking to the
    slipped intermediate):

        dG(L) = dG0 - alpha * (L - 40)      =>   phi = 1 / (1 + exp(dG/RT))

    Reparameterised with beta = alpha/RT (per repeat) and L50 = 40 + dG0/alpha
    (the length at which the slipped state is half-populated).
    """
    return 1.0 / (1.0 + np.exp(-beta * (L - L50)))


def step_L(L, step0=None, Ls=np.inf):
    """Repeats gained per expansion event at tract length L.

    The insertion copied during repair synthesis is the extrahelical slip-out
    itself (a (CTG)4 insertion was recovered from a 4-unit loop in vitro), so the
    step size is the slip-out size. Longer tracts support larger slip-outs, which
    is modelled as a linear growth with a characteristic scale Ls:

        s(L) = step0 * (1 + (L - 40) / Ls)

    Ls -> infinity recovers a length-independent step of `step0`.
    """
    s0 = LIT["step_repeats"] if step0 is None else step0
    return s0 * (1.0 + np.maximum(L - 40.0, 0.0) / Ls)


def clustered_fraction(L, L_sw=None, w_sw=None):
    """Fraction of slip-outs that occur in *clusters* rather than in isolation.

    Isolated short CTG/CAG slip-outs are repaired efficiently by hMutSbeta,
    whereas clustered slip-outs are repaired poorly (Panigrahi et al., PNAS
    2010), so the escape-to-expansion probability per slip-out rises sharply
    once a tract is long enough to host several slip-outs simultaneously.
    Modelled as a sigmoid switch in tract length.
    """
    if L_sw is None:
        raise ValueError("clustered_fraction: L_sw must be supplied; it is a "
                         "fitted parameter and has no default")
    w_sw = ASSUMED["w_switch"] if w_sw is None else w_sw
    return 1.0 / (1.0 + np.exp(-(np.asarray(L, float) - L_sw) / w_sw))


def lambda_L(L, A, beta, L50=None, msh3=1.0, fan1=1.0, footprint=None,
             R=None, L_sw=None, w_sw=None):
    """Per-cell expansion-event rate (events per year) at tract length L.

        lambda(L) = A * n(L) * phi(L) * [ (1 - C(L)) + R * C(L) ] * theta_MSH3 / fan1

    n(L)   MutSbeta site multiplicity (17-repeat footprint)
    phi(L) equilibrium occupancy of the extrahelical 1-unit slip-out state
    C(L)   fraction of slip-outs in poorly repaired clusters
    R      fold increase in expansion probability per clustered slip-out
           relative to an isolated one (>= 1)

    MutSbeta engagement is first-order in MSH3 abundance: one MutSbeta engaging
    one 1-unit slip-out is the productive species in >90% of CAG particles, so
    theta is taken proportional to relative MSH3 protein level. `fan1` is
    relative FAN1 activity; FAN1 opposes MutLgamma-dependent expansion, so the
    expansion flux is divided by it.
    """
    L50 = ASSUMED["L50_ceiling"] if L50 is None else L50
    if R is None:
        raise ValueError("lambda_L: R must be supplied; it is a fitted "
                         "parameter and has no default")
    C = clustered_fraction(L, L_sw, w_sw)
    return (A * n_sites(L, footprint) * phi_slipout(L, beta, L50)
            * ((1.0 - C) + R * C) * msh3 / fan1)


def lam(L, par, msh3=1.0, fan1=1.0):
    """Convenience wrapper: expansion rate at L for a fitted parameter dict."""
    return lambda_L(L, par["A"], par["beta"], par.get("L50"), msh3=msh3,
                    fan1=fan1, R=par.get("R"), L_sw=par.get("L_sw"),
                    w_sw=par.get("w_sw"))


def mean_transit_time(L1, L2, par, msh3=1.0, fan1=1.0, step0=None):
    """Deterministic (mean-field) time to expand from L1 to L2.

    Integrates dL/dt = s(L) * lambda(L), i.e.
        T = int_{L1}^{L2} dL / (s(L) * lambda(L)).
    Used for fast parameter fitting; verified against the exact Gillespie
    simulation in `03_coupling.py`.
    """
    grid = np.linspace(L1, L2, 4001)
    integrand = 1.0 / (step_L(grid, step0, par.get("Ls", np.inf))
                       * lam(grid, par, msh3=msh3, fan1=fan1))
    return float(np.trapezoid(integrand, grid))


def mean_first_passage(L0, L_target, par, n_cells=4000, seed=0, msh3=1.0,
                       fan1=1.0, step=None, contraction_ratio=None,
                       contraction_step=None, t_max=400.0):
    """Mean first-passage time from L0 to L_target by exact Gillespie.

    This is the stochastic counterpart of `mean_transit_time`, which integrates
    the mean-field trajectory dL/dt = s(L)*lambda(L). For a nonlinear lambda
    the two are not the same quantity, and contractions enter the stochastic
    process but not the integral, so the comparison is a real check rather
    than an identity. Elimination is not simulated: the elimination hazard is
    zero below the 150-repeat threshold, so it cannot affect any of the
    published transit-time anchors.

    Returns (mean, standard error, fraction of cells that reached L_target).
    """
    rng = np.random.default_rng(seed)
    step0 = LIT["step_repeats"] if step is None else step
    Ls_par = par.get("Ls", np.inf)
    cr = ASSUMED["contraction_ratio"] if contraction_ratio is None else contraction_ratio
    cs = ASSUMED["contraction_step"] if contraction_step is None else contraction_step
    times = []
    for _ in range(n_cells):
        t, L = 0.0, float(L0)
        while t < t_max and L < L_target:
            lam_exp = float(lam(L, par, msh3=msh3, fan1=fan1))
            tot = lam_exp * (1.0 + cr)
            t += rng.exponential(1.0 / tot)
            if rng.random() * tot < lam_exp:
                s_mean = max(float(step_L(L, step0, Ls_par)), 1.0 + 1e-9)
                L = min(L + float(rng.geometric(1.0 / s_mean)),
                        ASSUMED["L_max"])
            else:
                L = max(L - cs, 5.0)
        if L >= L_target:
            times.append(t)
    arr = np.asarray(times, float)
    if arr.size == 0:
        return float("nan"), float("nan"), 0.0
    return (float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(arr.size)),
            arr.size / n_cells)


def _params_from_vector(v):
    """(log lambda40, beta, log R, L_sw) -> full parameter dict, A derived."""
    loglam40, beta, logR, L_sw = v
    lam40 = float(np.exp(loglam40))
    R = float(np.exp(logR))
    L50 = ASSUMED["L50_ceiling"]
    w_sw = ASSUMED["w_switch"]
    C40 = float(clustered_fraction(40.0, L_sw, w_sw))
    A = lam40 / (n_sites(40.0) * phi_slipout(40.0, beta, L50)
                 * ((1 - C40) + R * C40))
    return {"A": float(A), "beta": float(beta), "L50": L50, "R": R,
            "L_sw": float(L_sw), "w_sw": w_sw, "lambda40": lam40,
            "Ls": np.inf}


def calibrate(targets=None, verbose=False, fold_target=None, assumed=None,
              quick=False):
    """Fit the rate law to four published SPN kinetic anchors.

    Fitted parameters (4): lambda(40), beta, R, L_switch.
    Targets (Handsaker et al., Cell 2025):
        T1  40 -> 60 CAG  in 50 years
        T2  60 -> 80 CAG  in 12 years
        T3  80 -> 150 CAG in ~4 years
        T4  lambda(90)/lambda(70) = 6   (reported as "at least 6-fold")

    T4 is reported as a lower bound; it is used here as an equality target, so
    the fitted model sits at the slowest acceleration consistent with the
    report. Remaining anchors (events per year below 55 CAG, fraction of life
    spent below 150 CAG, post-threshold survival, fraction of SPNs exceeding the
    inherited length) are held out for validation.
    """
    from scipy.optimize import least_squares

    if assumed:
        saved = {k: ASSUMED[k] for k in assumed}
        ASSUMED.update(assumed)
    else:
        saved = None

    t = dict(t1=LIT["t_40_60_years"], t2=LIT["t_60_80_years"],
             t3=LIT["t_80_150_years"],
             t4=LIT["rate_rise_70_90_min"] if fold_target is None else fold_target)
    if targets:
        t.update(targets)

    def resid(v):
        p = _params_from_vector(v)
        m = [mean_transit_time(40, 60, p),
             mean_transit_time(60, 80, p),
             mean_transit_time(80, 150, p),
             float(lam(90.0, p) / lam(70.0, p))]
        obs = [t["t1"], t["t2"], t["t3"], t["t4"]]
        return [np.log(mi / oi) for mi, oi in zip(m, obs)]

    starts = ([(-2.67, 0.006, np.log(17.7), 80.2)] if quick else
              [(a, b, c, d) for a in (-3.0, -2.0) for b in (0.02, 0.05, 0.10)
               for c in (np.log(3.0), np.log(10.0), np.log(50.0))
               for d in (70.0, 80.0, 90.0)])

    best = None
    for x0 in starts:
        try:
            r = least_squares(resid, list(x0),
                              bounds=([-12, 1e-4, 0.0, 50.0],
                                      [5, 1.0, np.log(1e5), 140.0]),
                              xtol=1e-15, ftol=1e-15, gtol=1e-15)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        raise RuntimeError("calibration failed for every starting point")

    par = _params_from_vector(best.x)
    par["cost"] = float(best.cost)
    if saved:
        ASSUMED.update(saved)
    if verbose:
        print("fitted:", {k: (round(v, 6) if isinstance(v, float) else v)
                          for k, v in par.items()})
    return par


# ---------------------------------------------------------------------------
# Layer 3 : PD link from MSH3 protein level to expansion suppression
# ---------------------------------------------------------------------------

def prevention_from_msh3(msh3_pct):
    """% prevention of somatic expansion as a function of striatal MSH3 protein
    (% of untreated control), from the published linear regression
    (Sci Rep 2024, doi:10.1038/s41598-024-52667-3: y-intercept 92.4, x-intercept 110.4, R^2 = 0.65).
    """
    y0 = LIT["pd_y_intercept"]
    x0 = LIT["pd_x_intercept"]
    return np.clip(y0 * (1.0 - np.asarray(msh3_pct, float) / x0), 0.0, 100.0)


def model_prevention_from_msh3(msh3_pct, order=1.0):
    """Model-predicted % prevention if the productive MutSbeta species requires
    `order` simultaneously bound molecules (rate ~ theta**order, theta ~ MSH3).
    order = 1 is the first-order model used in this work; order = 3 is the
    cooperative ("tandem-requirement") alternative tested in the model
    comparison.
    """
    m = np.asarray(msh3_pct, float) / 100.0
    return 100.0 * (1.0 - m ** order)


# ---------------------------------------------------------------------------
# Layer 2 : single-cell stochastic simulation
# ---------------------------------------------------------------------------

def death_hazard(L, L_death_50=None, width=None, h_max=None):
    """Elimination hazard (per year) for an SPN whose tract has reached L.

    Zero below the Phase-C threshold (150 CAG), then a sigmoid in L centred on
    L_death_50, so that the observed "months rather than decades" post-threshold
    decline emerges from the very fast expansion at long lengths rather than
    being imposed.
    """
    L_death_50 = ASSUMED["L_death_50"] if L_death_50 is None else L_death_50
    width = ASSUMED["death_width"] if width is None else width
    h_max = ASSUMED["h_death_max"] if h_max is None else h_max
    L = np.asarray(L, float)
    h = h_max / (1.0 + np.exp(-(L - L_death_50) / width))
    return np.where(L < LIT["L_phaseC"], 0.0, h)


def simulate_cell(L0, par, rng, t_max=110.0, msh3=1.0, fan1=1.0,
                  step=None, contraction_ratio=None, contraction_step=None,
                  hazard_kw=None, record_threshold=True,
                  t_treat=np.inf, msh3_treat=1.0):
    """Exact Gillespie simulation of one SPN.

    Events: expansion (+ geometric step), contraction (-contraction_step),
    elimination. Treatment is a step change in relative MSH3 protein from
    `msh3` to `msh3_treat` at age `t_treat`; the simulation advances exactly to
    that time and resamples, so the piecewise-constant rates stay exact.

    Returns dict with time of crossing 150 CAG (t_C), time of death (t_death),
    final and maximal tract length, and total years lived below 150 CAG.
    """
    step0 = LIT["step_repeats"] if step is None else step
    Ls_par = par.get("Ls", np.inf)
    cr = ASSUMED["contraction_ratio"] if contraction_ratio is None else contraction_ratio
    cs = ASSUMED["contraction_step"] if contraction_step is None else contraction_step
    hazard_kw = hazard_kw or {}

    t = 0.0
    L = float(L0)
    L_max_seen = L
    t_C = np.nan
    years_below = 0.0

    while t < t_max:
        m_now = msh3 if t < t_treat else msh3_treat
        lam_exp = float(lam(L, par, msh3=m_now, fan1=fan1))
        lam_con = cr * lam_exp
        h = float(death_hazard(L, **hazard_kw))
        tot = lam_exp + lam_con + h
        # Reached only at complete MSH3 ablation (msh3 = 0) below the
        # elimination threshold, where no event of any kind can occur. Skip
        # forward to the next change in treatment status instead of drawing
        # from an exponential with zero rate.
        if tot <= 0:
            t_next = min(t_treat if t < t_treat else t_max, t_max)
            if L < LIT["L_phaseC"]:
                years_below += t_next - t
            t = t_next
            if t >= t_max:
                break
            continue
        dt = rng.exponential(1.0 / tot)
        # do not step past a change in treatment status
        horizon = t_treat if t < t_treat else np.inf
        if t + dt > horizon:
            if L < LIT["L_phaseC"]:
                years_below += min(horizon, t_max) - t
            t = horizon
            continue
        if L < LIT["L_phaseC"]:
            years_below += min(dt, t_max - t)
        t += dt
        if t >= t_max:
            break
        u = rng.random() * tot
        if u < lam_exp:
            # repeats gained = size of the slip-out copied during repair synthesis;
            # drawn from a geometric distribution with mean s(L) (minimum 1)
            s_mean = max(float(step_L(L, step0, Ls_par)), 1.0 + 1e-9)
            gained = float(rng.geometric(1.0 / s_mean))
            L = min(L + gained, ASSUMED["L_max"])
            L_max_seen = max(L_max_seen, L)
            if record_threshold and np.isnan(t_C) and L >= LIT["L_phaseC"]:
                t_C = t
        elif u < lam_exp + lam_con:
            L = max(L - cs, 5.0)
        else:
            return {"t_C": t_C, "t_death": t, "L_final": L,
                    "L_max": L_max_seen, "years_below_150": years_below,
                    "died": True}

    return {"t_C": t_C, "t_death": np.nan, "L_final": L, "L_max": L_max_seen,
            "years_below_150": years_below, "died": False}


_CELL_DTYPE = np.dtype([("t_C", "f8"), ("t_death", "f8"), ("L_final", "f8"),
                        ("L_max", "f8"), ("years_below_150", "f8"), ("died", "?")])


def simulate_population(L0, par, n_cells=2000, seed=0, **kw):
    """Simulate a population of independent SPNs; returns a structured array."""
    rng = np.random.default_rng(seed)
    arr = np.empty(n_cells, dtype=_CELL_DTYPE)
    for i in range(n_cells):
        o = simulate_cell(L0, par, rng, **kw)
        arr[i] = (o["t_C"], o["t_death"], o["L_final"], o["L_max"],
                  o["years_below_150"], o["died"])
    return arr


def simulate_population_tau(L0, par, n_cells=20000, dt=1.0 / 24.0, t_max=110.0,
                            msh3=1.0, fan1=1.0, t_treat=np.inf, msh3_treat=1.0,
                            k_nca=0.0, contraction_ratio=None,
                            contraction_step=None, step=None, hazard_kw=None,
                            seed=0, return_curve=True):
    """Vectorised tau-leaping simulation of a coupled SPN population.

    Identical to `simulate_population` in the limit dt -> 0 and k_nca = 0, but
    fast enough for parameter sweeps and able to carry a *non-cell-autonomous*
    term: the elimination hazard of every surviving SPN is multiplied by
    (1 + k_nca * f_lost), where f_lost is the fraction of the SPN population
    already eliminated. k_nca = 0 recovers the purely cell-autonomous model.

    Expansion within a step is exact in the step-size distribution: the number
    of events is Poisson(lambda*dt) and the total repeats gained is the sum of
    that many geometric draws, i.e. a negative-binomial variate.
    """
    rng = np.random.default_rng(seed)
    step0 = LIT["step_repeats"] if step is None else step
    Ls_par = par.get("Ls", np.inf)
    cr = ASSUMED["contraction_ratio"] if contraction_ratio is None else contraction_ratio
    cs = ASSUMED["contraction_step"] if contraction_step is None else contraction_step
    hazard_kw = hazard_kw or {}

    L = np.full(n_cells, float(L0))
    alive = np.ones(n_cells, dtype=bool)
    t_death = np.full(n_cells, np.nan)
    t_C = np.full(n_cells, np.nan)
    L_max = L.copy()
    years_below = np.zeros(n_cells)

    n_steps = int(np.ceil(t_max / dt))
    times = np.empty(n_steps + 1)
    surviving = np.empty(n_steps + 1)
    times[0] = 0.0
    surviving[0] = 1.0

    for k in range(n_steps):
        t = k * dt
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            times[k + 1:] = t + dt * np.arange(1, n_steps - k + 1)
            surviving[k + 1:] = 0.0
            break
        m_now = msh3 if t < t_treat else msh3_treat
        Li = L[idx]
        years_below[idx] += dt * (Li < LIT["L_phaseC"])

        lamv = lam(Li, par, msh3=m_now, fan1=fan1)
        n_ev = rng.poisson(lamv * dt)
        n_con = rng.poisson(cr * lamv * dt)

        gain = np.zeros(idx.size)
        pos = n_ev > 0
        if np.any(pos):
            s_mean = np.maximum(step_L(Li[pos], step0, Ls_par), 1.0 + 1e-9)
            p = 1.0 / s_mean
            gain[pos] = n_ev[pos] + rng.negative_binomial(n_ev[pos], p)
        Lnew = np.clip(Li + gain - n_con * cs, 5.0, ASSUMED["L_max"])

        newly = np.isnan(t_C[idx]) & (Lnew >= LIT["L_phaseC"])
        t_C[idx[newly]] = t + dt
        L[idx] = Lnew
        L_max[idx] = np.maximum(L_max[idx], Lnew)

        f_lost = 1.0 - alive.mean()
        h = death_hazard(Lnew, **hazard_kw) * (1.0 + k_nca * f_lost)
        p_die = 1.0 - np.exp(-h * dt)
        dies = rng.random(idx.size) < p_die
        if np.any(dies):
            d = idx[dies]
            alive[d] = False
            t_death[d] = t + dt

        times[k + 1] = t + dt
        surviving[k + 1] = alive.mean()

    arr = np.empty(n_cells, dtype=_CELL_DTYPE)
    arr["t_C"] = t_C
    arr["t_death"] = t_death
    arr["L_final"] = L
    arr["L_max"] = L_max
    arr["years_below_150"] = years_below
    arr["died"] = ~alive
    if return_curve:
        return arr, times, surviving
    return arr


def survival_curve(arr, ages):
    """Fraction of SPNs still alive at each age."""
    d = arr["t_death"]
    d = d[~np.isnan(d)]
    n = len(arr)
    return np.array([1.0 - np.sum(d <= a) / n for a in ages])


def age_at_loss(arr, frac, ages=None):
    """Age at which a given fraction of SPNs has been lost (NaN if never)."""
    ages = np.linspace(0, 110, 2201) if ages is None else ages
    s = survival_curve(arr, ages)
    tgt = 1.0 - frac
    idx = np.where(s <= tgt)[0]
    return float(ages[idx[0]]) if len(idx) else np.nan


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def cap_score(age, cag):
    return age * (cag - LIT["cap_offset"])


def ranen_predicted_onset(cag, parent_aao=None):
    """External sanity anchor: Ranen 'analogous model' regression reported by
    Orth & Schwenke (PLOS Curr Huntingt Dis 2011)."""
    p = LIT["ranen_mean_parent_aao"] if parent_aao is None else parent_aao
    return (LIT["ranen_intercept"] + LIT["ranen_parent_coef"] * p
            + LIT["ranen_cag_coef"] * cag)
