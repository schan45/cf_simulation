"""Editable research priors and example inputs; none are patient-validated.

Tuple tables use (example, lower tuning bound, upper tuning bound).
Bounds are exploratory, NOT confidence intervals. See modell_modszertan.md.
"""

import numpy as np

from mikrobiom_modell import (
    DailyDose, Ecology, Model, Patient, PD, PK, immune_clearance,
)

SPECIES = {
    "PA": "Pseudomonas aeruginosa",
    "BCC": "Burkholderia cepacia complex",
    "SA": "Staphylococcus aureus",
    "HI": "Haemophilus influenzae",
    "SM": "Stenotrophomonas maltophilia",
    "VD": "Veillonella dispar",
    "PM": "Prevotella melaninogenica",
    "SS": "Streptococcus sanguinis",
    "FN": "Fusobacterium nucleatum",
    "RM": "Rothia mucilaginosa",
    "AX": "Achromobacter xylosoxidans",
}

# Task-1 priors retained for comparability, not reinterpreted as fitted values.
R_TABLE = {
    "PA": (1.3, .05, 1.66), "BCC": (.5, .02, 1), "SA": (.66, .014, 1.02),
    "HI": (.5, .02, 1.4), "SM": (.99, .005, 1.2), "VD": (.5, .01, .8),
    "PM": (.25, .01, .7), "SS": (.42, .02, .6), "FN": (.2, .005, .3),
    "RM": (.3, .01, .8), "AX": (.42, .02, 1),
}
K_TABLE = {
    "PA": (5e8, 1e7, 1e10), "BCC": (1e8, 1e6, 1e9), "SA": (5e7, 1e5, 1e9),
    "HI": (1e7, 1e5, 1e8), "SM": (1e6, 1e4, 1e8), "VD": (1e5, 1e3, 1e7),
    "PM": (3e6, 1e4, 1e8), "SS": (1e8, 1e6, 1e9), "FN": (1e3, 1e1, 1e5),
    "RM": (1e7, 1e5, 1e8), "AX": (1e7, 1e5, 1e8),
}

# (source j, recipient i, a_ij, lower, upper); source leads in mikrobiom_parameterek.md.
# Many edges are association/analogy hypotheses; use INTERACTION_SCALE=0 as a control.
INTERACTIONS = [
    ("PA", "SA", -.55, -.90, -.30), ("SA", "PA", .20, .05, .40),
    ("PA", "BCC", -.65, -.90, -.50), ("BCC", "PA", .15, -.30, .40),
    ("PA", "HI", -.15, -.30, -.05), ("HI", "PA", .10, .00, .25),
    ("PA", "SM", -.65, -.90, -.50), ("SM", "PA", .20, .05, .40),
    ("PA", "AX", -.20, -.30, .00), ("AX", "PA", -.55, -.75, -.35),
    ("PA", "VD", -.25, -.30, -.05), ("VD", "PA", .10, .00, .20),
    ("PA", "PM", -.40, -.50, -.10), ("PM", "PA", .25, .10, .45),
    ("PA", "SS", -.10, -.10, .40), ("SS", "PA", -.50, -.90, -.20),
    ("PA", "FN", -.10, -.25, -.05), ("FN", "PA", .30, .10, .50),
    ("PA", "RM", -.40, -.50, -.10), ("RM", "PA", .20, .05, .40),
    ("SA", "BCC", .15, .00, .40), ("SM", "BCC", .10, .00, .30),
    ("SA", "HI", .15, -.30, .30), ("HI", "SA", -.05, -.30, .00),
    ("SA", "SM", -.05, -.15, .00), ("AX", "SM", -.25, -.50, -.10),
    ("SS", "SA", -.30, -.50, -.20), ("FN", "SA", .35, .10, .50),
    ("SA", "FN", .30, .10, .50), ("VD", "SA", -.05, -.30, .00),
    ("SA", "VD", -.15, -.30, .00), ("PM", "SA", -.15, -.30, .00),
    ("SA", "PM", -.15, -.30, .00), ("RM", "SA", -.05, -.30, .00),
    ("SS", "HI", -.15, -.30, -.05), ("SS", "VD", .30, .30, .80),
    ("VD", "SS", .20, .20, .60), ("FN", "SS", .10, .10, .60),
    ("SS", "FN", -.45, -.50, .20), ("VD", "FN", .35, .30, .80),
    ("FN", "VD", .25, .20, .60), ("PM", "FN", .25, .10, .50),
    ("FN", "PM", .20, .10, .50), ("VD", "PM", .10, .05, .30),
    ("PM", "VD", .10, .05, .30), ("PM", "SS", .15, .10, .35),
    ("SS", "PM", .15, .10, .35), ("SS", "RM", .20, .10, .30),
    ("RM", "SS", .10, .10, .30), ("FN", "RM", .05, .05, .30),
    ("RM", "FN", .05, .05, .30), ("VD", "RM", .05, .05, .30),
    ("RM", "VD", .05, .05, .30), ("PM", "RM", .10, .05, .30),
    ("RM", "PM", .05, .05, .30), ("VD", "BCC", .10, .00, .25),
    ("PM", "BCC", .10, .00, .25), ("SS", "BCC", .25, .05, .35),
    ("FN", "BCC", .15, .05, .30), ("BCC", "VD", -.20, -.25, .00),
    ("BCC", "PM", -.40, -.50, -.10), ("BCC", "SS", -.35, -.35, -.05),
    ("BCC", "FN", -.15, -.25, .00), ("VD", "SM", .10, .00, .25),
    ("PM", "SM", .10, .00, .25), ("SS", "SM", .30, .10, .40),
    ("FN", "SM", .15, -.05, .30), ("SM", "VD", -.15, -.25, .00),
    ("SM", "PM", -.15, -.30, .05), ("SM", "SS", .05, .00, .30),
    ("SM", "FN", -.20, -.30, .00), ("SM", "RM", .10, .00, .30),
    ("VD", "AX", .05, .00, .25), ("PM", "AX", .05, .00, .25),
    ("SS", "AX", .10, .05, .35), ("FN", "AX", .10, .05, .30),
]
INTERACTION_SCALE = 1.0

PK_PARAMETERS = {
    "tobramycin": PK("inhaled", 2, .10, 1.0, 5.2,
                    "Geller 2002 DOI 10.1378/chest.122.1.219: sputum scale; "
                    "local half-life and transit are unvalidated priors."),
    "aztreonam": PK("inhaled", 2, .10, 1.0, 12.2,
                   "EMA Cayston section 5.2: 726 ug/g at 10 min after 75 mg; "
                   "1 g/mL sputum density assumed; local half-life/transit are priors."),
    "azithromycin": PK("oral", 2, 2.0, 72.0, .01,
                      "Wilms 2006 DOI 10.1097/01.ftd.0000195617.69721.a5: "
                      "12-53 mg/L sputum at chronic 500 mg/day; all local PK values are priors."),
}
PK_TUNING = {
    "tobramycin": {"mean_transit_h": (.02, .5), "local_half_life_h": (.25, 4),
                  "dose_to_local_per_l": (.5, 20)},
    "aztreonam": {"mean_transit_h": (.02, .5), "local_half_life_h": (.25, 4),
                 "dose_to_local_per_l": (1, 40)},
    "azithromycin": {"mean_transit_h": (.5, 6), "local_half_life_h": (24, 168),
                    "dose_to_local_per_l": (.001, .05)},
}

# Per row: kill_max [1/h], inhibition_max [0,1], EC50 [mg/L].
# Zeros mean explicitly unmodelled action, NOT proven biological neutrality.
PD_ROWS = {
    "PA":  {"tobramycin": (2.0, 0, 100), "aztreonam": (1.8, 0, 100), "azithromycin": (0, 0, 1)},
    "BCC": {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
    "SA":  {"tobramycin": (1.0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, .85, 1)},
    "HI":  {"tobramycin": (1.0, 0, 100), "aztreonam": (1.0, 0, 100), "azithromycin": (0, .85, 1)},
    "SM":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
    "VD":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
    "PM":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
    "SS":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, .80, 1)},
    "FN":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
    "RM":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, .80, 1)},
    "AX":  {"tobramycin": (0, 0, 100), "aztreonam": (0, 0, 100), "azithromycin": (0, 0, 1)},
}
PD_TUNING = {
    "kill_max_h": (0, 5), "inhibition_max": (0, 1),
    "ec50_mg_l": (.01, 1000), "hill": (1, 4),
    "resistance_fold": (1, 1000), "initial_resistant_fraction": (0, .1),
    "mutation_per_division": (0, 1e-5), "fitness_cost": (0, .3),
}
PD_PARAMETERS = {
    species: {
        drug: PD(
            kill_max_h=row[0], inhibition_max=row[1], ec50_mg_l=row[2],
            hill=2, resistance_fold=16,
            initial_resistant_fraction=1e-4 if row[0] or row[1] else 0,
            mutation_per_division=1e-8 if row[0] or row[1] else 0,
            fitness_cost=.05 if row[0] or row[1] else 0,
            provenance="Unfitted demonstration prior; see species/drug evidence table in modell_modszertan.md.",
        )
        for drug, row in rows.items()
    }
    for species, rows in PD_ROWS.items()
}

IMMUNE_BASE_H = .01
LUNG_WEIGHTS = {"mild": 1.0, "moderate": .7, "severe": .4}
DISEASE_WEIGHTS = {
    "pancreatic_insufficiency": .05, "cf_related_diabetes": .15,
    "impaired_glucose_tolerance": .08, "liver_disease": .05,
    "bone_disease": .02, "sinus_disease": .02, "recurrent_dios": .03,
}
# ALL immune weights are user-requested hypothetical scenarios, not clinical estimates.
PATIENT = Patient(lung_category="moderate", comorbidities=("pancreatic_insufficiency",))

DURATION_DAYS = 14.0
OUTPUT_STEP_H = .5
MAX_STEP_H = .5
SNAPSHOT_DAYS = (0, 3, 7, 10, 14)
DETECTION_LIMIT_CFU_ML = 1.0

SELECTED_SPECIES = ("PA", "SA", "SS")
INITIAL_MODE = "relative_anchor"  # alternatives: "log10" or "cfu"
INITIAL_PERCENTAGES = {"PA": 65.0, "SA": 20.0, "SS": 13.0}
ANCHOR_SPECIES = "PA"
ANCHOR_LOG10_CFU_ML = float(np.log10(6.5e7))
INITIAL_LOG10_CFU_ML = {"PA": float(np.log10(6.5e7)), "SA": float(np.log10(2e7)),
                        "SS": float(np.log10(1.3e7))}
INITIAL_CFU_ML = {"PA": 6.5e7, "SA": 2e7, "SS": 1.3e7}
OTHERS_CFU_ML = 2e6  # Used only in log10/cfu input modes.

# Each row represents one calendar day. Modify/replace rows to change individual days.
# No treatment is inferred from missing rows. Input is daily TOTAL mg, not mg per dose.
REGIMEN = (
    tuple(DailyDose("aztreonam", day, 750, 3) for day in range(1, 15))
    + tuple(DailyDose("azithromycin", day, 300, 2) for day in range(8, 15))
)
# Tobramycin is implemented but OFF in the example. Example syntax for a single day:
# DailyDose("tobramycin", day=1, daily_mg=600, administrations=2, hours=(8, 20))


def build_model(selected: tuple[str, ...] = SELECTED_SPECIES, patient: Patient = PATIENT) -> Model:
    if not selected or len(set(selected)) != len(selected) or set(selected) - set(SPECIES):
        raise ValueError("Select distinct species from SPECIES.")
    if not np.isfinite(INTERACTION_SCALE) or INTERACTION_SCALE < 0:
        raise ValueError("INTERACTION_SCALE must be finite and nonnegative.")
    index = {s: i for i, s in enumerate(selected)}
    a = -np.eye(len(selected))
    seen = set()
    for source, recipient, value, lower, upper in INTERACTIONS:
        if source not in SPECIES or recipient not in SPECIES or source == recipient:
            raise ValueError("Invalid directed interaction.")
        if (source, recipient) in seen or not lower <= value <= upper:
            raise ValueError("Duplicate interaction or example outside tuning bounds.")
        seen.add((source, recipient))
        if source in index and recipient in index:
            a[index[recipient], index[source]] = value * INTERACTION_SCALE
    for table in (R_TABLE, K_TABLE):
        for s in selected:
            if not 0 < table[s][1] <= table[s][0] <= table[s][2]:
                raise ValueError(f"{s}: parameter outside positive tuning interval.")
    ecology = Ecology(
        selected, np.array([R_TABLE[s][0] for s in selected]),
        np.array([K_TABLE[s][0] for s in selected]), a,
    )
    clearance = immune_clearance(patient, IMMUNE_BASE_H, LUNG_WEIGHTS, DISEASE_WEIGHTS)
    return Model(ecology, PK_PARAMETERS, {s: PD_PARAMETERS[s] for s in selected}, clearance)
