"""A magyar felület adatellenőrzése és eredménytáblái."""

import math
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from mikrobiom_modell import DRUGS, DailyDose, Result, initial_from_relative

NAME = "Baktérium"
PERCENT = "Relatív abundancia (%)"
CFU = "CFU/mL"
LOG = "log10 CFU/mL"
SELECTED = "Kiválasztva"
DRUG = "Antibiotikum"
DOSE = "Napi összadag (mg)"
COUNT = "Beadások / nap"

DRUG_LABELS = {
    "tobramycin": "Tobramycin (inhalált)",
    "aztreonam": "Aztreonam (inhalált)",
    "azithromycin": "Azithromycin (orális)",
}
LUNG_LABELS = {
    "mild": "Enyhe – FEV1 ≥70%",
    "moderate": "Közepes – 40% ≤ FEV1 <70%",
    "severe": "Súlyos – FEV1 <40%",
}
DISEASE_LABELS = {
    "pancreatic_insufficiency": "Hasnyálmirigy-elégtelenség",
    "cf_related_diabetes": "CF-hez kapcsolódó diabetes",
    "impaired_glucose_tolerance": "Csökkent glükóztolerancia",
    "liver_disease": "Májbetegség",
    "bone_disease": "Csontbetegség",
    "sinus_disease": "Orrmelléküreg-betegség",
    "recurrent_dios": "Visszatérő disztális bélelzáródás (DIOS)",
}


def _number(value: object, label: str) -> float:
    if not isinstance(value, (int, float, np.number)) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label}: szám megadása szükséges.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label}: véges szám megadása szükséges.")
    return number


def microbiome_input(
    table: pd.DataFrame, selected: Sequence[str], *, relative: bool, log_scale: bool,
    others_absolute: float = 0,
) -> tuple[dict[str, float], float]:
    if not selected or len(set(selected)) != len(selected) or len(table) != len(selected):
        raise ValueError("Válassz legalább egy különböző baktériumot.")
    column = LOG if log_scale else CFU
    absolute: dict[str, float] = {}
    for species, (_, row) in zip(selected, table.iterrows()):
        raw = row[column]
        if pd.isna(raw):
            continue
        value = _number(raw, f"{species} – {column}")
        if log_scale:
            if not -300 <= value <= 300:
                raise ValueError("A log10 CFU/mL értéke −300 és 300 között legyen.")
            value = 10 ** value
        elif value < 0:
            raise ValueError("A CFU/mL nem lehet negatív.")
        absolute[species] = value
    if relative:
        percentages = {
            species: _number(value, f"{species} – százalék")
            for species, value in zip(selected, table[PERCENT])
        }
        if any(v < 0 for v in percentages.values()) or math.fsum(percentages.values()) > 100:
            raise ValueError("A százalékok legyenek nemnegatívak, összegük legfeljebb 100%.")
        anchors = [s for s in selected if absolute.get(s, 0) > 0 and percentages[s] > 0]
        if not anchors:
            raise ValueError("Legalább egy pozitív százalékú baktériumhoz pozitív abszolút CFU/mL szükséges.")
        anchor = anchors[0]
        initial, others = initial_from_relative(percentages, anchor, math.log10(absolute[anchor]))
        for species, value in absolute.items():
            if not math.isclose(value, initial[species], rel_tol=0.01, abs_tol=0):
                raise ValueError(
                    f"{species}: az abszolút mennyiség és a százalék ellentmond egymásnak "
                    "(1%-os kerekítési tűrés). Hagyj meg egy horgonyt, vagy javítsd az adatokat."
                )
    else:
        if set(absolute) != set(selected):
            raise ValueError("Abszolút módban minden kiválasztott baktérium mennyiségét add meg.")
        initial, others = absolute, _number(others_absolute, "Others CFU/mL")
        if others < 0:
            raise ValueError("Az Others mennyisége nem lehet negatív.")
    if not math.isfinite(math.fsum(initial.values()) + others) or math.fsum(initial.values()) + others <= 0:
        raise ValueError("A teljes kezdeti baktériumterhelés legyen pozitív és véges.")
    return initial, others


def regimen_input(table: pd.DataFrame, days: int) -> tuple[DailyDose, ...]:
    if len(table) != len(DRUGS):
        raise ValueError("A kezelési táblázatnak három gyógyszersort kell tartalmaznia.")
    rows = []
    for drug, (_, row) in zip(DRUGS, table.iterrows()):
        if not row[SELECTED]:
            continue
        active = [d for d in range(1, days + 1) if row[f"{d}. nap"]]
        if not active:
            raise ValueError(f"{DRUG_LABELS[drug]}: jelölj ki legalább egy kezelési napot.")
        mg = _number(row[DOSE], DOSE)
        count = _number(row[COUNT], COUNT)
        if mg <= 0 or count < 1 or not count.is_integer():
            raise ValueError(f"{DRUG_LABELS[drug]}: pozitív mg/nap és pozitív egész beadásszám szükséges.")
        rows.extend(DailyDose(drug, d, mg, int(count)) for d in active)
    return tuple(rows)


def snapshot_input(text: str, duration: float) -> tuple[float, ...]:
    try:
        values = tuple(sorted(set(float(s.strip().replace(",", ".")) for s in text.split(";"))))
    except ValueError as exc:
        raise ValueError("Az időpontokat pontosvesszővel válaszd el, például: 0; 2,5; 7; 14.") from exc
    if not values or any(not math.isfinite(v) or not 0 <= v <= duration for v in values):
        raise ValueError("Minden időpontnak a 0. nap és a szimuláció vége közé kell esnie.")
    return values


def final_table(result: Result, names: Mapping[str, str], control: Result | None = None) -> pd.DataFrame:
    start = np.append(result.populations[0], result.others_cfu_ml)
    end = np.append(result.populations[-1], result.others_cfu_ml)
    change_pct = np.divide(100 * (end - start), start, out=np.full_like(start, np.nan), where=start > 0)
    delta_log = np.full_like(start, np.nan)
    valid = (start > 0) & (end > 0)
    delta_log[valid] = np.log10(end[valid]) - np.log10(start[valid])
    table = pd.DataFrame({
        NAME: [names[s] for s in result.species] + ["Others (állandó háttér)"],
        "Kezdeti CFU/mL": start, "Végső CFU/mL": end,
        "Változás (CFU/mL)": end - start, "Változás (%)": change_pct,
        "Δlog10 (nyers)": delta_log, "Végső arány (%)": result.percentages[-1],
    })
    if control is not None:
        if control.species != result.species or control.time_h[-1] != result.time_h[-1]:
            raise ValueError("A kontroll taxonjai és időhorizontja nem egyeznek.")
        control_end = np.append(control.populations[-1], control.others_cfu_ml)
        table["Kontroll végső CFU/mL"] = control_end
        table["Kezelt − kontroll (CFU/mL)"] = end - control_end
    return table
