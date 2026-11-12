"""UI input tests and an end-to-end Streamlit session."""

from pathlib import Path
import unittest
from io import BytesIO
import zipfile

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from streamlit.testing.v1 import AppTest

import modell_parameterek as p
from mikrobiom_modell import dose_events, simulate
from ui_adatok import (
    CFU, COUNT, DOSE, DRUG, LOG, NAME, PERCENT, SELECTED,
    final_table, microbiome_input, regimen_input, snapshot_input,
)


class InputTests(unittest.TestCase):
    def relative_table(self) -> pd.DataFrame:
        return pd.DataFrame({
            NAME: ["PA", "SA", "SS"], CFU: [6.5e7, np.nan, np.nan], PERCENT: [65., 20., 13.],
        })

    def test_relative_anchor_and_others(self) -> None:
        initial, others = microbiome_input(self.relative_table(), ("PA", "SA", "SS"), relative=True, log_scale=False)
        assert_allclose(list(initial.values()), [6.5e7, 2e7, 1.3e7])
        self.assertAlmostEqual(others, 2e6)

    def test_log_mode_and_multiple_anchors(self) -> None:
        table = self.relative_table().rename(columns={CFU: LOG})
        table[LOG] = [np.log10(6.5e7), np.log10(2e7), np.nan]
        initial, _ = microbiome_input(table, ("PA", "SA", "SS"), relative=True, log_scale=True)
        self.assertAlmostEqual(initial["SS"], 1.3e7)
        table.loc[1, LOG] = 4
        with self.assertRaisesRegex(ValueError, "ellentmond"):
            microbiome_input(table, ("PA", "SA", "SS"), relative=True, log_scale=True)

    def test_invalid_relative_input(self) -> None:
        for percentages, absolute in (([70, 20, 13], 6.5e7), ([65, 20, 13], np.nan),
                                      ([0, 20, 13], 6.5e7), ([-1, 20, 13], 6.5e7)):
            table = self.relative_table()
            table[PERCENT] = percentages
            table.loc[0, CFU] = absolute
            with self.assertRaises(ValueError):
                microbiome_input(table, ("PA", "SA", "SS"), relative=True, log_scale=False)

    def test_absolute_mode(self) -> None:
        table = self.relative_table()
        table[CFU] = [1e6, 0, 2e6]
        initial, others = microbiome_input(table, ("PA", "SA", "SS"), relative=False,
                                          log_scale=False, others_absolute=3e6)
        self.assertEqual(initial["SA"], 0)
        self.assertEqual(others, 3e6)
        table.loc[1, CFU] = np.nan
        with self.assertRaises(ValueError):
            microbiome_input(table, ("PA", "SA", "SS"), relative=False, log_scale=False)

    def test_regimen(self) -> None:
        table = pd.DataFrame([
            {SELECTED: enabled, DRUG: drug, DOSE: mg, COUNT: count,
             **{f"{d}. nap": (drug == "aztreonam" or drug == "azithromycin" and d >= 8)
                for d in range(1, 15)}}
            for enabled, drug, mg, count in (
                (False, "tobramycin", 600, 2), (True, "aztreonam", 750, 3), (True, "azithromycin", 300, 2),
            )
        ])
        events = dose_events(regimen_input(table, 14), 14)
        self.assertEqual(len(events), 56)
        self.assertEqual(sum(e.mg for e in events), 750 * 14 + 300 * 7)
        table[SELECTED] = False
        self.assertEqual(regimen_input(table, 14), ())
        table.loc[0, SELECTED] = True
        with self.assertRaisesRegex(ValueError, "napot"):
            regimen_input(table, 14)
        table.loc[0, "1. nap"] = True
        table.loc[0, COUNT] = 0
        with self.assertRaisesRegex(ValueError, "beadásszám"):
            regimen_input(table, 14)

    def test_snapshots_and_final_cfu(self) -> None:
        self.assertEqual(snapshot_input("0; 2,5; 14; 2.5", 14), (0., 2.5, 14.))
        for text in ("", "nan", "15", "1, 2"):
            with self.assertRaises(ValueError):
                snapshot_input(text, 14)
        result = simulate(p.build_model(("PA",)), {"PA": 1e6}, (), duration_days=.1, others_cfu_ml=1e5)
        table = final_table(result, p.SPECIES, result)
        assert_allclose(table["Végső CFU/mL"], [result.populations[-1, 0], 1e5])
        assert_allclose(table["Kezelt − kontroll (CFU/mL)"], 0)


class AppTests(unittest.TestCase):
    def test_full_default_run_and_stale_results(self) -> None:
        app = AppTest.from_file(str(Path(__file__).with_name("app.py")), default_timeout=120).run()
        self.assertFalse(app.exception, str(app.exception))
        self.assertFalse(app.error, str(app.error))
        app.button(key="run").click().run()
        self.assertFalse(app.exception, str(app.exception))
        self.assertFalse(app.error, str(app.error))
        completed = app.session_state["completed"]
        self.assertEqual(completed["days"], 14)
        self.assertEqual(len(completed["summary"]), 4)
        with zipfile.ZipFile(BytesIO(completed["zip"])) as archive:
            for name in ("vegeredmeny.csv", "kontroll_CFU_mL.csv", "futas_beallitasai.json",
                         "szimulacios_gorbek.png", "mikrobiom_idopontok.png"):
                self.assertIn(name, archive.namelist())
        app.number_input(key="days").set_value(10).run()
        self.assertFalse(app.exception)
        self.assertTrue(app.error)
        self.assertTrue(app.button(key="run").disabled)
        self.assertFalse(app.tabs)
        app.text_input(key="snapshots").set_value("0; 3; 7").run()
        self.assertFalse(app.error)
        self.assertFalse(app.tabs)

    def test_invalid_profile_and_empty_taxa(self) -> None:
        app = AppTest.from_file(str(Path(__file__).with_name("app.py")), default_timeout=30).run()
        app.multiselect(key="species").set_value([]).run()
        self.assertFalse(app.exception)
        self.assertTrue(app.error)
        self.assertTrue(app.button(key="run").disabled)


if __name__ == "__main__":
    unittest.main()
