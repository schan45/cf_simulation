"""Fixed equations and input validation. Units: h, CFU/mL, mg, mg/L.

Research model, not a clinical dose calculator. See modell_modszertan.md.
"""

from dataclasses import dataclass
from itertools import product
from collections.abc import Mapping, Sequence
import math
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.special import expit

Array = NDArray[np.float64]
DRUGS = ("tobramycin", "aztreonam", "azithromycin")
PHENOTYPES = np.array(list(product((0, 1), repeat=3)), dtype=int)
COMORBIDITIES = (
    "pancreatic_insufficiency", "cf_related_diabetes",
    "impaired_glucose_tolerance", "liver_disease", "bone_disease",
    "sinus_disease", "recurrent_dios",
)


def _finite(value: float, name: str, *, positive: bool = False) -> None:
    if not math.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}.")


@dataclass(frozen=True)
class Patient:
    sex: str = "female"
    age_years: float = 22
    bmi_kg_m2: float = 21
    lung_category: str = "moderate"
    comorbidities: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.sex not in ("female", "male"):
            raise ValueError("sex must be female or male.")
        _finite(self.age_years, "age_years")
        _finite(self.bmi_kg_m2, "bmi_kg_m2", positive=True)
        if self.lung_category not in ("mild", "moderate", "severe"):
            raise ValueError("Unknown FEV1 category.")
        if len(set(self.comorbidities)) != len(self.comorbidities):
            raise ValueError("Duplicate comorbidities.")
        if set(self.comorbidities) - set(COMORBIDITIES):
            raise ValueError("Unknown comorbidity.")
        if {"cf_related_diabetes", "impaired_glucose_tolerance"} <= set(self.comorbidities):
            raise ValueError("Select diabetes or impaired glucose tolerance, not both.")


def immune_clearance(
    patient: Patient, base_h: float, lung_weights: Mapping[str, float],
    disease_weights: Mapping[str, float],
) -> float:
    _finite(base_h, "immune base")
    if set(lung_weights) != {"mild", "moderate", "severe"}:
        raise ValueError("Provide all three lung weights.")
    if set(disease_weights) != set(COMORBIDITIES):
        raise ValueError("Provide all comorbidity weights.")
    for name, value in (*lung_weights.items(), *disease_weights.items()):
        _finite(value, name)
    return base_h * lung_weights[patient.lung_category] * math.exp(
        -sum(disease_weights[d] for d in patient.comorbidities)
    )


@dataclass(frozen=True)
class DailyDose:
    drug: str
    day: int
    daily_mg: float
    administrations: int
    hours: tuple[float, ...] | None = None

    def __post_init__(self) -> None:
        if self.drug not in DRUGS:
            raise ValueError(f"Unknown drug: {self.drug}")
        if isinstance(self.day, bool) or not isinstance(self.day, int) or self.day < 1:
            raise ValueError("day must be a positive integer (day 1 starts at t=0).")
        _finite(self.daily_mg, "daily_mg")
        if isinstance(self.administrations, bool) or not isinstance(self.administrations, int):
            raise ValueError("administrations must be an integer.")
        if self.daily_mg == 0:
            if self.administrations != 0 or self.hours not in (None, ()):
                raise ValueError("A treatment-free day requires 0 mg and 0 administrations.")
        elif self.administrations < 1:
            raise ValueError("Positive daily dose needs at least one administration.")
        if self.hours is not None and self.daily_mg > 0:
            if len(self.hours) != self.administrations or len(set(self.hours)) != len(self.hours):
                raise ValueError("hours must contain one distinct time per administration.")
            if any(not math.isfinite(h) or not 0 <= h < 24 for h in self.hours):
                raise ValueError("Administration hours must be in [0,24).")


@dataclass(frozen=True)
class DoseEvent:
    drug: str
    time_h: float
    mg: float


def dose_events(regimen: Sequence[DailyDose], duration_days: float) -> tuple[DoseEvent, ...]:
    _finite(duration_days, "duration_days", positive=True)
    seen: set[tuple[str, int]] = set()
    events = []
    for row in regimen:
        key = (row.drug, row.day)
        if key in seen:
            raise ValueError(f"Duplicate drug/day row: {key}")
        seen.add(key)
        if row.day - 1 >= duration_days:
            raise ValueError(f"Day {row.day} starts outside the simulation.")
        if row.daily_mg == 0:
            continue
        hours = row.hours
        if hours is None:
            hours = tuple(24 * k / row.administrations for k in range(row.administrations))
        for hour in hours:
            t = 24 * (row.day - 1) + hour
            if t >= duration_days * 24:
                raise ValueError("A dose falls outside the simulation; shorten the regimen.")
            events.append(DoseEvent(row.drug, t, row.daily_mg / row.administrations))
    return tuple(sorted(events, key=lambda e: (e.time_h, e.drug)))


@dataclass(frozen=True)
class PK:
    route: str
    transit_stages: int
    mean_transit_h: float
    local_half_life_h: float
    dose_to_local_per_l: float
    provenance: str

    def __post_init__(self) -> None:
        if self.route not in ("inhaled", "oral"):
            raise ValueError("Unknown route.")
        if (isinstance(self.transit_stages, bool)
                or not isinstance(self.transit_stages, int) or self.transit_stages < 2):
            raise ValueError("Gamma/Erlang transit_stages must be an integer >= 2.")
        for name in ("mean_transit_h", "local_half_life_h", "dose_to_local_per_l"):
            _finite(getattr(self, name), name, positive=True)


@dataclass(frozen=True)
class PD:
    kill_max_h: float
    inhibition_max: float
    ec50_mg_l: float
    hill: float
    resistance_fold: float
    initial_resistant_fraction: float
    mutation_per_division: float
    fitness_cost: float
    provenance: str

    def __post_init__(self) -> None:
        _finite(self.kill_max_h, "kill_max_h")
        for name in ("ec50_mg_l", "hill", "resistance_fold"):
            _finite(getattr(self, name), name, positive=True)
        if self.resistance_fold < 1:
            raise ValueError("Resistance must not reduce EC50.")
        for name in ("inhibition_max", "initial_resistant_fraction",
                     "mutation_per_division", "fitness_cost"):
            value = getattr(self, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"{name} must be in [0,1].")


@dataclass(frozen=True)
class Ecology:
    species: tuple[str, ...]
    r_h: Array
    capacity_cfu_ml: Array
    interaction: Array

    def __post_init__(self) -> None:
        n = len(self.species)
        if n == 0 or len(set(self.species)) != n:
            raise ValueError("Select at least one distinct species.")
        for value, shape in (
            (self.r_h, (n,)), (self.capacity_cfu_ml, (n,)), (self.interaction, (n, n)),
        ):
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ValueError(f"Invalid ecology array, expected shape {shape}.")
        if np.any(self.r_h <= 0) or np.any(self.capacity_cfu_ml <= 0):
            raise ValueError("r and K must be positive.")
        if not np.all(np.diag(self.interaction) == -1):
            raise ValueError("A[i,i] must equal -1.")


def initial_from_relative(
    percentages: Mapping[str, float], anchor_species: str, anchor_log10_cfu_ml: float,
) -> tuple[dict[str, float], float]:
    if not percentages or "Others" in percentages:
        raise ValueError("Provide selected taxa only; Others is calculated as the remainder.")
    for key, value in percentages.items():
        _finite(value, f"percentage {key}")
    total = math.fsum(percentages.values())
    if total > 100:
        raise ValueError("Selected percentages must sum to at most 100.")
    if anchor_species not in percentages or percentages[anchor_species] <= 0:
        raise ValueError("The absolute anchor must have a positive relative abundance.")
    if not math.isfinite(anchor_log10_cfu_ml) or not -300 <= anchor_log10_cfu_ml <= 300:
        raise ValueError("Invalid log10 CFU/mL anchor.")
    total_cfu = 10 ** anchor_log10_cfu_ml / (percentages[anchor_species] / 100)
    if not math.isfinite(total_cfu):
        raise ValueError("Initial total CFU/mL overflow.")
    return {k: total_cfu * v / 100 for k, v in percentages.items()}, total_cfu * (100 - total) / 100


def initial_from_log10(log_values: Mapping[str, float]) -> dict[str, float]:
    if not log_values:
        raise ValueError("Provide at least one taxon.")
    if any(not math.isfinite(v) or not -300 <= v <= 300 for v in log_values.values()):
        raise ValueError("Invalid log10 CFU/mL. Use direct CFU/mL input for exact zero.")
    return {k: 10 ** v for k, v in log_values.items()}


def relative_abundance(populations: Array, others_cfu_ml: float) -> Array:
    _finite(others_cfu_ml, "Others")
    if populations.ndim != 2 or np.any(populations < 0) or not np.all(np.isfinite(populations)):
        raise ValueError("Populations must be a finite nonnegative time-by-species array.")
    values = np.column_stack((populations, np.full(len(populations), others_cfu_ml)))
    totals = values.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("Relative abundance is undefined for zero total biomass.")
    return values / totals * 100


def log10_observation(populations: Array, detection_limit: float) -> Array:
    _finite(detection_limit, "detection_limit", positive=True)
    if np.any(populations < 0) or not np.all(np.isfinite(populations)):
        raise ValueError("Cannot observe invalid populations.")
    return np.log10(np.maximum(populations, detection_limit))


@dataclass(frozen=True)
class Model:
    ecology: Ecology
    pk: Mapping[str, PK]
    pd: Mapping[str, Mapping[str, PD]]
    immune_h: float

    def __post_init__(self) -> None:
        _finite(self.immune_h, "immune_h")
        if set(self.pk) != set(DRUGS):
            raise ValueError("Specify PK for exactly the three supported antibiotics.")
        if self.pk["azithromycin"].route != "oral" or any(
            self.pk[d].route != "inhaled" for d in DRUGS[:2]
        ):
            raise ValueError("Tobramycin/aztreonam must be inhaled; azithromycin must be oral.")
        if set(self.pd) != set(self.ecology.species):
            raise ValueError("PD species must match the selected ecology.")
        for species in self.ecology.species:
            if set(self.pd[species]) != set(DRUGS):
                raise ValueError(f"Explicit PD, including null assumptions, required for {species}.")
            if sum(self.pd[species][d].mutation_per_division for d in DRUGS) > 1:
                raise ValueError("The sum of mutation probabilities per division must not exceed 1.")


@dataclass(frozen=True)
class Result:
    time_h: Array
    species: tuple[str, ...]
    phenotype_cfu_ml: Array
    concentration_mg_l: Array
    others_cfu_ml: float
    events: tuple[DoseEvent, ...]
    roundoff_corrections: int

    @property
    def populations(self) -> Array:
        return self.phenotype_cfu_ml.sum(axis=2)

    @property
    def percentages(self) -> Array:
        return relative_abundance(self.populations, self.others_cfu_ml)

    @property
    def resistant_fractions(self) -> Array:
        numerator = self.phenotype_cfu_ml @ PHENOTYPES
        total = self.populations[:, :, None]
        return np.divide(numerator, total, out=np.full_like(numerator, np.nan), where=total > 0)


class _System:
    def __init__(self, model: Model) -> None:
        self.model = model
        self.n = len(model.ecology.species)
        self.pk_slices: dict[str, slice] = {}
        offset = 8 * self.n
        for drug in DRUGS:
            length = model.pk[drug].transit_stages + 1
            self.pk_slices[drug] = slice(offset, offset + length)
            offset += length
        self.size = offset
        self.cost = np.ones((self.n, 8))
        self.ec50 = np.zeros((self.n, 8, 3))
        self.kill = np.zeros((self.n, 1, 3))
        self.inhibit = np.zeros((self.n, 1, 3))
        self.hill = np.zeros((self.n, 1, 3))
        self.mutation = np.zeros((self.n, 3))
        self.initial_p = np.ones((self.n, 8))
        for i, species in enumerate(model.ecology.species):
            for j, drug in enumerate(DRUGS):
                pd = model.pd[species][drug]
                bit = PHENOTYPES[:, j]
                self.cost[i] *= (1 - pd.fitness_cost) ** bit
                self.ec50[i, :, j] = pd.ec50_mg_l * pd.resistance_fold ** bit
                self.kill[i, 0, j] = pd.kill_max_h
                self.inhibit[i, 0, j] = pd.inhibition_max
                self.hill[i, 0, j] = pd.hill
                self.mutation[i, j] = pd.mutation_per_division
                f = pd.initial_resistant_fraction
                self.initial_p[i] *= np.where(bit == 1, f, 1 - f)
        off = model.ecology.interaction.copy()
        np.fill_diagonal(off, 0)
        self.facilitation = np.maximum(off, 0)
        self.competition = -np.minimum(off, 0)

    def rhs(self, _time_h: float, state: Array) -> Array:
        out = np.zeros_like(state)
        concentrations = np.empty(3)
        for j, drug in enumerate(DRUGS):
            p = self.model.pk[drug]
            block = state[self.pk_slices[drug]]
            dblock = out[self.pk_slices[drug]]
            rate = p.transit_stages / p.mean_transit_h
            dblock[0] = -rate * block[0]
            dblock[1:-1] = rate * (block[:-2] - block[1:-1])
            dblock[-1] = p.dose_to_local_per_l * rate * block[-2] - math.log(2) / p.local_half_life_h * block[-1]
            concentrations[j] = block[-1]
        # Solvers may evaluate just outside the positive cone at internal stages.
        c = np.maximum(concentrations, 0)[None, None, :]
        with np.errstate(divide="ignore"):
            log_ratio = self.hill * (np.log(c) - np.log(self.ec50))
        effect = expit(log_ratio)
        static = np.prod(1 - self.inhibit * effect, axis=2)
        death = np.sum(self.kill * effect, axis=2)
        y = state[:8 * self.n].reshape(self.n, 8)
        ecology = self.model.ecology
        normalized = y.sum(axis=1) / ecology.capacity_cfu_ml
        birth = ecology.r_h[:, None] * (1 + self.facilitation @ normalized)[:, None]
        birth = birth * self.cost * static
        loss = ecology.r_h * (normalized + self.competition @ normalized) + self.model.immune_h
        dy = (birth - loss[:, None] - death) * y
        # Mutations occur during positive birth events, independently of exposure.
        for j in range(3):
            source = np.flatnonzero(PHENOTYPES[:, j] == 0)
            target = source + 2 ** (2 - j)
            flux = self.mutation[:, j, None] * birth[:, source] * y[:, source]
            dy[:, source] -= flux
            dy[:, target] += flux
        out[:8 * self.n] = dy.ravel()
        return out


def simulate(
    model: Model, initial_cfu_ml: Mapping[str, float], regimen: Sequence[DailyDose],
    *, duration_days: float = 14, output_step_h: float = 0.5, max_step_h: float = 0.5,
    snapshot_days: Sequence[float] = (), others_cfu_ml: float = 0,
    rtol: float = 1e-7, atol: float = 1e-9,
) -> Result:
    for value, name in ((output_step_h, "output_step_h"), (max_step_h, "max_step_h"),
                        (rtol, "rtol"), (atol, "atol")):
        _finite(value, name, positive=True)
    _finite(others_cfu_ml, "others_cfu_ml")
    events = dose_events(regimen, duration_days)
    if set(initial_cfu_ml) != set(model.ecology.species):
        raise ValueError("Initial taxa must exactly match selected species.")
    for key, value in initial_cfu_ml.items():
        _finite(value, key)
    if math.fsum(initial_cfu_ml.values()) + others_cfu_ml <= 0:
        raise ValueError("Initial total biomass must be positive.")
    if any(not math.isfinite(t) or not 0 <= t <= duration_days for t in snapshot_days):
        raise ValueError("Snapshot times must be within the simulation.")
    end = 24 * duration_days
    rise_times = [
        e.time_h + model.pk[e.drug].mean_transit_h * factor
        for e in events for factor in (.5, 1, 2, 4)
        if e.time_h + model.pk[e.drug].mean_transit_h * factor < end
    ]
    times = np.unique(np.concatenate((
        np.arange(0, end, output_step_h), [end],
        np.array(snapshot_days) * 24, [e.time_h for e in events], rise_times,
    )))
    system = _System(model)
    state = np.zeros(system.size)
    initial = np.array([initial_cfu_ml[s] for s in model.ecology.species])
    state[:8 * system.n] = (initial[:, None] * system.initial_p).ravel()
    samples = np.empty((len(times), system.size))
    samples[0] = state
    by_time: dict[float, list[DoseEvent]] = {}
    for e in events:
        by_time.setdefault(e.time_h, []).append(e)
    bounds = sorted({0.0, end, *by_time})
    corrections = 0
    for left, right in zip(bounds[:-1], bounds[1:]):
        for event in by_time.get(left, []):
            state[system.pk_slices[event.drug].start] += event.mg
        selected = np.flatnonzero((times > left) & (times <= right))
        solution = solve_ivp(
            system.rhs, (left, right), state, method="LSODA",
            t_eval=times[selected], max_step=max_step_h, rtol=rtol, atol=atol,
        )
        if not solution.success:
            raise RuntimeError(f"ODE integration failed at {left} h: {solution.message}")
        block = solution.y.T
        if not np.all(np.isfinite(block)) or np.any(block < -50 * atol):
            raise FloatingPointError("Nonfinite or materially negative state; review parameters/tolerances.")
        corrections += int(np.count_nonzero(block < 0))
        block = np.maximum(block, 0)
        samples[selected] = block
        state = block[-1].copy()
    if corrections:
        warnings.warn(
            f"{corrections} negative solver-roundoff entries (<=50*atol) set to zero; "
            "tighten atol for low-count inference.", RuntimeWarning, stacklevel=2,
        )
    concentrations = np.column_stack([samples[:, system.pk_slices[d].stop - 1] for d in DRUGS])
    result = Result(
        times, model.ecology.species, samples[:, :8 * system.n].reshape(-1, system.n, 8),
        concentrations, others_cfu_ml, events, corrections,
    )
    # Also checks the observation denominator, without hiding undefined composition.
    result.percentages
    return result
