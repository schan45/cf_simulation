"""Önálló, magyar nyelvű Streamlit felület a TDK kutatási modellhez."""

from dataclasses import asdict
from io import BytesIO
import hashlib
import importlib
import json
import logging
from pathlib import Path
import tempfile
import warnings
import zipfile

import numpy as np
import pandas as pd
import streamlit as st

import modell_parameterek as p
from mikrobiom_modell import DRUGS, Patient, relative_abundance, simulate
from szimulacio_futtatas import save_results
from ui_adatok import (
    CFU, COUNT, DISEASE_LABELS, DOSE, DRUG, DRUG_LABELS, LOG, LUNG_LABELS,
    NAME, PERCENT, SELECTED, final_table, microbiome_input, regimen_input, snapshot_input,
)

importlib.reload(p)
ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title="CF-mikrobiom | TDK", page_icon="🧬", layout="wide")
st.title("CF-mikrobiom szimulátor")
st.caption("TDK kutatási felület")

with st.sidebar:
    st.header("Futtatási beállítások")
    days = int(st.number_input("Szimuláció hossza (nap)", min_value=1, value=14, step=1, key="days"))
    step = st.selectbox("Kimeneti lépés (óra)", (.1, .25, .5, 1.0), index=2)
    snapshots_text = st.text_input(
        "Összetétel vizsgálata (napok)", value="0; 3; 7; 10; 14",
        help="Pontosvesszővel elválasztva; tört nap is megadható.", key="snapshots",
    )
    with_control = st.checkbox("Gyógyszermentes kontroll is", value=True, key="control")
    st.caption("A kontroll ugyanazt a mikrobiomot és immunprofilt használja, antibiotikum nélkül.")

st.header("1. Páciensprofil")
col1, col2, col3, col4 = st.columns(4)
sex = col1.selectbox("Nem", ("female", "male"), format_func=lambda v: {"female": "Nő", "male": "Férfi"}[v])
age = col2.number_input("Életkor (év)", min_value=0, value=22, step=1)
bmi = col3.number_input("BMI (kg/m²)", min_value=.1, value=21.0, step=.1)
lung = col4.selectbox("Tüdőbetegség kategóriája", tuple(LUNG_LABELS), index=1, format_func=LUNG_LABELS.get)
diseases = st.multiselect("CF-hez kapcsolódó betegségek", tuple(DISEASE_LABELS),
                         default=["pancreatic_insufficiency"], format_func=DISEASE_LABELS.get)
st.caption("A nem, életkor és BMI csak metaadat. A FEV1-kategória és a társbetegségek "
           "demonstrációs, nem validált súlyokkal módosítják az immunkonstanst.")

st.header("2. Mikrobiomprofil")
selected = tuple(st.multiselect("Baktériumok kiválasztása", tuple(p.SPECIES),
                               default=list(p.SELECTED_SPECIES), format_func=p.SPECIES.get, key="species"))
left, right = st.columns(2)
mode = left.radio("Beviteli mód", ("Relatív abundancia + abszolút horgony", "Minden mennyiség abszolút"),
                  key="mode")
scale = right.radio("Abszolút értékek skálája", (CFU, LOG), horizontal=True, key="scale")
relative = mode.startswith("Relatív")
log_scale = scale == LOG
st.caption("Relatív módban minden kiválasztott faj százalékát és legalább egy pozitív abszolút "
           "mennyiséget adj meg. A többi abszolút cella üresen hagyható. Több horgonynak egyeznie kell "
           "(1% kerekítési tűrés). A becslés nem egyenértékű tenyésztéses méréssel.")
defaults = []
for species in selected:
    count = p.INITIAL_CFU_ML.get(species, p.K_TABLE[species][0] * .1)
    value = float(np.log10(count)) if log_scale else float(count)
    defaults.append({
        NAME: p.SPECIES[species], scale: value if not relative or species == "PA" else None,
        PERCENT: p.INITIAL_PERCENTAGES.get(species, 0.0),
    })
micro_df = pd.DataFrame(defaults, columns=[NAME, scale, PERCENT])
micro_df[scale] = micro_df[scale].astype(float)
micro_df[PERCENT] = micro_df[PERCENT].astype(float)
editor_key = f"micro_{','.join(selected)}_{relative}_{log_scale}"
micro = st.data_editor(
    micro_df, hide_index=True, width="stretch", key=editor_key,
    disabled=[NAME] + ([PERCENT] if not relative else []),
    column_config={
        scale: st.column_config.NumberColumn(scale, min_value=-300.0 if log_scale else 0.0,
                                             format="%.5f" if log_scale else "%.4e"),
        PERCENT: st.column_config.NumberColumn(PERCENT, min_value=0, max_value=100, format="%.3f"),
    },
)
others_absolute = 0.0
if not relative:
    others_absolute = st.number_input("Others – abszolút háttér (CFU/mL)", min_value=0.0, value=2e6, format="%.4e")
st.caption("Others = a 100%-hoz hiányzó arány. Az ebből kapott abszolút háttér a futás alatt állandó; "
           "aránya a többi taxon változásával módosul. Az ábrákon szürke.")

validation_errors = []
try:
    initial, others = microbiome_input(
        micro, selected, relative=relative, log_scale=log_scale, others_absolute=others_absolute,
    )
except ValueError as exc:
    validation_errors.append(str(exc))
else:
    preview = pd.DataFrame({
        NAME: [p.SPECIES[s] for s in selected] + ["Others"],
        "Számított CFU/mL": list(initial.values()) + [others],
        "Számított arány (%)": relative_abundance(np.array([list(initial.values())]), others)[0],
    })
    st.dataframe(preview.style.format({"Számított CFU/mL": "{:.5e}", "Számított arány (%)": "{:.3f}"}))

st.header("3. Antibiotikumprofil")
st.caption("A napi mg az ÖSSZES napi adag: mg/beadás = mg/nap ÷ beadásszám. "
           "A pipák jelölik a kezelési napokat. A beadások egyenletesen oszlanak el a napban, "
           "az első a nap kezdetén. Az 1. nap t=0, a 8. nap t=7 nap időpontban kezdődik.")
regimen_df = pd.DataFrame([
    {
        SELECTED: d != "tobramycin", DRUG: DRUG_LABELS[d],
        DOSE: {"tobramycin": 600.0, "aztreonam": 750.0, "azithromycin": 300.0}[d],
        COUNT: {"tobramycin": 2, "aztreonam": 3, "azithromycin": 2}[d],
        **{f"{day}. nap": (d == "aztreonam" or d == "azithromycin" and day >= 8)
           for day in range(1, days + 1)},
    } for d in DRUGS
])
regimen_table = st.data_editor(
    regimen_df, hide_index=True, width="stretch", disabled=[DRUG], key=f"regimen_{days}",
    column_config={
        SELECTED: st.column_config.CheckboxColumn(SELECTED),
        DOSE: st.column_config.NumberColumn(DOSE, min_value=0.0),
        COUNT: st.column_config.NumberColumn(COUNT, min_value=0, step=1),
        **{f"{d}. nap": st.column_config.CheckboxColumn(str(d), width="small") for d in range(1, days + 1)},
    },
)
st.caption("Időhorizont vagy mikrobiom-beviteli mód váltásakor az érintett táblázat alapértékei "
           "újratöltődnek. Ellenőrizd a pipákat és a kezdeti mennyiségeket futtatás előtt.")
if "cf_related_diabetes" in diseases and "impaired_glucose_tolerance" in diseases:
    validation_errors.append("Diabetes és csökkent glükóztolerancia közül itt csak egyet válassz.")
try:
    regimen = regimen_input(regimen_table, days)
    snapshot_days = snapshot_input(snapshots_text, days)
except ValueError as exc:
    validation_errors.append(str(exc))

model = None
if not validation_errors:
    try:
        patient = Patient(sex, age, bmi, lung, tuple(diseases))
        model = p.build_model(selected, patient)
    except ValueError as exc:
        validation_errors.append(f"A modellbeállítások hibásak: {exc}")
    else:
        st.metric("Számított immunkonstans", f"{model.immune_h:.6f} h⁻¹")

for error in validation_errors:
    st.error(error)

signature = None
if not validation_errors and model is not None:
    signature = hashlib.sha256(json.dumps({
        "model": asdict(model), "patient": asdict(patient), "initial": initial, "others": others,
        "regimen": [asdict(r) for r in regimen], "days": days, "snapshots": snapshot_days,
        "step": step, "control": with_control, "detection": p.DETECTION_LIMIT_CFU_ML,
    }, default=lambda a: a.tolist(), sort_keys=True).encode()).hexdigest()

if st.button("Szimuláció futtatása", type="primary", disabled=bool(validation_errors), key="run"):
    st.session_state.pop("completed", None)
    if model is not None:
        try:
            with st.spinner("Szimuláció és ábrák készítése…"), warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                solver = {"output_step_h": step, "max_step_h": .5, "rtol": 1e-7, "atol": 1e-9}
                result = simulate(model, initial, regimen, duration_days=days,
                                  snapshot_days=snapshot_days, others_cfu_ml=others, **solver)
                control = simulate(model, initial, (), duration_days=days,
                                   others_cfu_ml=others, **solver) if with_control else None
                summary = final_table(result, p.SPECIES, control)
                with tempfile.TemporaryDirectory(prefix="tdk_ui_") as folder:
                    output = Path(folder)
                    save_results(
                        result, model, patient, initial, regimen, snapshot_days, p.DETECTION_LIMIT_CFU_ML,
                        output, solver, "UI: becsült CFU/mL relatív horgonnyal" if relative else "UI: felhasználói abszolút bemenet",
                    )
                    summary.to_csv(output / "vegeredmeny.csv", index=False, encoding="utf-8-sig")
                    if control is not None:
                        control_data = pd.DataFrame(control.populations, columns=control.species)
                        control_data.insert(0, "time_days", control.time_h / 24)
                        control_data["Others_CFU_mL"] = others
                        control_data.to_csv(output / "kontroll_CFU_mL.csv", index=False)
                    files = {path.name: path.read_bytes() for path in output.iterdir() if path.is_file()}
                buffer = BytesIO()
                with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                    for name, content in files.items():
                        archive.writestr(name, content)
                st.session_state.completed = {
                    "signature": signature, "files": files, "zip": buffer.getvalue(), "summary": summary,
                    "days": days, "warnings": list(dict.fromkeys(str(w.message) for w in caught)),
                }
        except (ValueError, RuntimeError, FloatingPointError, OverflowError, OSError) as exc:
            logging.exception("A szimuláció vagy az export sikertelen")
            st.error(f"A futás sikertelen. Ellenőrizd a bemeneteket és a modellparamétereket. Részlet: {exc}")

completed = st.session_state.get("completed")
if completed is not None and completed["signature"] != signature:
    st.info("A beállítások megváltoztak. A korábbi eredményt elrejtettük; futtasd újra a szimulációt.")
elif completed is not None:
    st.header("4. Eredmények")
    st.success(f"A {completed['days']} napos kutatási szimuláció elkészült.")
    if completed["warnings"]:
        with st.expander("Numerikus figyelmeztetések", expanded=True):
            st.warning("A solver kerekítési maradékokra figyelmeztethet. Ez nem biológiai kihalás. "
                       "A részletek az alábbiakban láthatók.")
            for warning in completed["warnings"]:
                st.code(warning, language=None)
    curves, composition, endpoint, downloads = st.tabs(
        ["Szimulációs görbék", "Összetétel és rezisztencia", "Végső abszolút mennyiségek", "Letöltés"]
    )
    with curves:
        st.image(completed["files"]["szimulacios_gorbek.png"])
    with composition:
        st.image(completed["files"]["mikrobiom_idopontok.png"])
        st.image(completed["files"]["rezisztens_frakciok.png"])
    with endpoint:
        summary = completed["summary"]
        total_start = summary["Kezdeti CFU/mL"].sum()
        total_end = summary["Végső CFU/mL"].sum()
        a, b = st.columns(2)
        a.metric("Kezdeti teljes terhelés (Others-szel)", f"{total_start:.4e} CFU/mL")
        b.metric("Végső teljes terhelés (Others-szel)", f"{total_end:.4e} CFU/mL",
                 delta=f"{100 * (total_end / total_start - 1):+.2f}%", delta_color="off")
        st.dataframe(summary, hide_index=True, width="stretch",
                     column_config={c: st.column_config.NumberColumn(c, format="%.5e")
                                    for c in summary if "CFU/mL" in c})
        st.caption("Az abszolút érték koncentráció (CFU/mL), nem a tüdő teljes CFU-ja. "
                   "A nem értelmezhető változások (például 0 kezdeti CFU-nál) üresek.")
        st.subheader("Kutatási értelmezés")
        st.write("A táblázat megmutatja, mely taxon mennyisége nőtt vagy csökkent, és – ha kértél "
                 "kontrollt – miben tér el a gyógyszermentes forgatókönyvtől.")
        oral = {"VD", "PM", "SS", "FN", "RM"}
        present = summary[summary[NAME].isin([p.SPECIES[s] for s in oral])]
        if not present.empty:
            st.write("Kiválasztott orális/kommensális csoportba sorolt taxonok:")
            st.dataframe(present[[NAME, "Végső CFU/mL", "Végső arány (%)"]], hide_index=True)
        st.info("Nincs általános „jó baktérium” címke: a hatás környezet-, faj- és törzsfüggő. "
                "A klinikai sikerhez betegadatok, tünetek, légzésfunkció és validált célértékek kellenek. ")
    with downloads:
        st.download_button("Teljes eredménycsomag (ZIP)", completed["zip"], "tdk_szimulacio.zip",
                           mime="application/zip")
        st.download_button("Végső CFU/mL táblázat (CSV)", completed["files"]["vegeredmeny.csv"],
                           "vegeredmeny.csv", mime="text/csv")
        st.caption("A ZIP tartalmazza az ábrákat, idősorokat, dózisokat, fenotípusállapotokat és a tényleges "
                   "futási beállításokat. A betegprofil is bekerül: csak azonosító nélküli adatot adj meg.")
