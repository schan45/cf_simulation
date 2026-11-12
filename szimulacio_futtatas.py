"""Run the editable example or call run_simulation with explicit model inputs."""

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
from collections.abc import Mapping, Sequence
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

import modell_parameterek as parameters
from mikrobiom_modell import (
    DRUGS, PHENOTYPES, DailyDose, Model, Patient, Result, initial_from_log10,
    initial_from_relative, log10_observation, simulate,
)

RESEARCH_WARNING = (
    "Research demonstration only: unfitted ecological/PD/immune priors and uncertain "
    "local PK. No clinical dose recommendation or validated patient prediction."
)


def _json_value(value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def load_regimen_csv(path: Path) -> tuple[DailyDose, ...]:
    """UTF-8 CSV columns: drug,day,daily_mg,administrations,hours (optional; separated by ;)."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not {
            "drug", "day", "daily_mg", "administrations"
        } <= set(reader.fieldnames):
            raise ValueError("Missing required regimen CSV columns.")
        rows = []
        for line, row in enumerate(reader, 2):
            try:
                hours_text = row.get("hours", "") or ""
                rows.append(DailyDose(
                    row["drug"], int(row["day"]), float(row["daily_mg"]),
                    int(row["administrations"]),
                    tuple(float(h) for h in hours_text.split(";")) if hours_text else None,
                ))
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Invalid regimen CSV row {line}: {exc}") from exc
    return tuple(rows)


def _snapshot_indices(result: Result, days: Sequence[float]) -> list[int]:
    if not days:
        raise ValueError("Select at least one snapshot.")
    indices = []
    for day in days:
        matches = np.flatnonzero(np.isclose(result.time_h, day * 24, rtol=0, atol=1e-9))
        if len(matches) != 1:
            raise ValueError(f"Snapshot {day} days was not sampled; include it in simulate().")
        indices.append(int(matches[0]))
    return indices


def plot_results(result: Result, snapshot_days: Sequence[float], detection_limit: float) -> dict[str, Figure]:
    days = result.time_h / 24
    colors = {s: plt.get_cmap("tab20")(list(parameters.SPECIES).index(s)) for s in result.species}
    drug_colors = {"tobramycin": "#0072B2", "aztreonam": "#D55E00", "azithromycin": "#CC79A7"}
    observed = log10_observation(result.populations, detection_limit)
    labels = [parameters.SPECIES[s] for s in result.species] + ["Others (allando CFU/mL)"]
    palette = [colors[s] for s in result.species] + ["#bbbbbb"]
    fig, axes = plt.subplots(4, 1, figsize=(13, 13), sharex=True, layout="constrained")
    for i, s in enumerate(result.species):
        axes[0].plot(days, observed[:, i], color=colors[s], label=parameters.SPECIES[s])
        # Censored observation change, not a claim about extinction.
        axes[2].plot(days, observed[:, i] - observed[0, i], color=colors[s])
    if observed.min() <= np.log10(detection_limit) + .5:
        axes[0].axhline(np.log10(detection_limit), color=".5", linestyle=":", linewidth=.8)
    axes[0].set_ylabel("log10 CFU/mL\n(detektálási küszöbbel)")
    axes[0].legend(fontsize=8, ncols=2, loc="upper right")
    axes[0].set_title("Szimulációs görbék", fontsize=13)
    active = {e.drug for e in result.events}
    for j, drug in enumerate(DRUGS):
        if drug in active:
            axes[1].plot(days, result.concentration_mg_l[:, j], color=drug_colors[drug], label=drug)
            event_days = [e.time_h / 24 for e in result.events if e.drug == drug]
            axes[3].plot(event_days, np.full(len(event_days), -.025 - .018 * j), "|",
                         transform=axes[3].get_xaxis_transform(), color=drug_colors[drug],
                         clip_on=False, markersize=5)
    if active:
        axes[1].set_yscale("symlog", linthresh=.1)
        axes[1].legend(fontsize=9)
    else:
        axes[1].text(.5, .5, "Nincs antibiotikum", transform=axes[1].transAxes, ha="center")
    axes[1].set_ylabel("Helyi C [mg/L]\n(symlog, kuszob 0.1)")
    axes[2].axhline(0, color=".5", linestyle=":", linewidth=.8)
    axes[2].set_ylabel("Változás a kezdethez\nDelta log10 CFU/mL")
    axes[3].stackplot(days, result.percentages.T, labels=labels, colors=palette)
    axes[3].set_ylim(0, 100)
    axes[3].set_ylabel("Relatív abundancia [%]")
    axes[3].legend(loc="upper center", ncols=2, fontsize=8)
    axes[3].set_xlabel("Idő [nap]; alsó színes jelek: dózisok")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.set_xlim(days[0], days[-1])

    indices = _snapshot_indices(result, snapshot_days)
    bars, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    bottom = np.zeros(len(indices))
    x = np.arange(len(indices))
    for i, (label, color) in enumerate(zip(labels, palette)):
        heights = result.percentages[indices, i]
        ax.bar(x, heights, bottom=bottom, color=color, label=label, width=.65)
        bottom += heights
    ax.set_xticks(x, [f"{d:g}. nap" for d in snapshot_days])
    ax.set_ylim(0, 100)
    ax.set_ylabel("Relatív abundancia [%]")
    ax.set_title("Mikrobiom összetétel a kiválasztott időpontokban")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=9)

    resistance, axes_r = plt.subplots(3, 1, figsize=(11, 8), sharex=True, layout="constrained")
    fractions = result.resistant_fractions
    for j, drug in enumerate(DRUGS):
        for i, species in enumerate(result.species):
            axes_r[j].plot(days, fractions[:, i, j] * 100, color=colors[species],
                           label=parameters.SPECIES[species])
        axes_r[j].set_ylabel(f"{drug}\nR frakcio [%]")
        axes_r[j].set_ylim(0, 100)
        axes_r[j].grid(alpha=.2)
    axes_r[0].set_title("Gyógyszer-specifikus rezisztens fenotipusok (nem MIC-meres)")
    axes_r[0].legend(fontsize=8, ncols=2)
    axes_r[-1].set_xlabel("Idő [nap]")
    return {"szimulacios_gorbek": fig, "mikrobiom_idopontok": bars, "rezisztens_frakciok": resistance}


def save_results(
    result: Result, model: Model, patient: Patient, initial_cfu_ml: Mapping[str, float],
    regimen: Sequence[DailyDose], snapshot_days: Sequence[float], detection_limit: float,
    output_dir: Path, solver_settings: Mapping[str, float], initial_provenance: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, fig in plot_results(result, snapshot_days, detection_limit).items():
        fig.savefig(output_dir / f"{name}.png", dpi=180)
        plt.close(fig)
    pct = result.percentages
    log_counts = log10_observation(result.populations, detection_limit)
    fractions = result.resistant_fractions
    with (output_dir / "idosor.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["time_h", "time_days"]
            + [f"{s}_CFU_mL" for s in result.species]
            + [f"{s}_log10_CFU_mL_censored" for s in result.species]
            + [f"{s}_percent" for s in (*result.species, "Others")]
            + [f"{d}_mg_L" for d in DRUGS]
            + [f"{s}_{d}_R_fraction" for s in result.species for d in DRUGS]
        )
        for k, t in enumerate(result.time_h):
            writer.writerow(
                [t, t / 24, *result.populations[k], *log_counts[k], *pct[k],
                 *result.concentration_mg_l[k],
                 *[v if np.isfinite(v) else "" for v in fractions[k].ravel()]]
            )
    with (output_dir / "idopontok.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_days", *result.species, "Others"])
        for day, index in zip(snapshot_days, _snapshot_indices(result, snapshot_days)):
            writer.writerow([day, *pct[index]])
    with (output_dir / "dozisok.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["drug", "time_h", "time_days", "dose_mg"])
        writer.writerows((e.drug, e.time_h, e.time_h / 24, e.mg) for e in result.events)
    np.savez_compressed(
        output_dir / "allapotok.npz", time_h=result.time_h,
        species=np.array(result.species), drugs=np.array(DRUGS), phenotype_bits=PHENOTYPES,
        phenotype_cfu_ml=result.phenotype_cfu_ml,
        concentration_mg_l=result.concentration_mg_l, others_cfu_ml=result.others_cfu_ml,
    )
    import scipy
    manifest = {
        "warning": RESEARCH_WARNING, "patient": asdict(patient), "model": asdict(model),
        "initial_cfu_ml": dict(initial_cfu_ml), "initial_provenance": initial_provenance,
        "others_cfu_ml": result.others_cfu_ml, "regimen": [asdict(row) for row in regimen],
        "snapshot_days": list(snapshot_days), "duration_days": float(result.time_h[-1] / 24),
        "detection_limit_cfu_ml": detection_limit, "solver": "LSODA segmented at dose events",
        "solver_settings": dict(solver_settings), "roundoff_corrections": result.roundoff_corrections,
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
    }
    (output_dir / "futas_beallitasai.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=_json_value, allow_nan=False),
        encoding="utf-8",
    )


def run_simulation(
    model: Model, initial_cfu_ml: Mapping[str, float], regimen: Sequence[DailyDose],
    *, patient: Patient, others_cfu_ml: float = 0, duration_days: float = 14,
    snapshot_days: Sequence[float] = (0, 7, 14), output_dir: Path = Path("eredmenyek"),
    output_step_h: float = .5, max_step_h: float = .5, detection_limit: float = 1,
    initial_provenance: str = "user-supplied; not asserted to be measured",
) -> Result:
    warnings.warn(RESEARCH_WARNING, UserWarning, stacklevel=2)
    if len(snapshot_days) == 0:
        raise ValueError("Select at least one snapshot.")
    if not np.isfinite(detection_limit) or detection_limit <= 0:
        raise ValueError("detection_limit must be positive.")
    solver = {"output_step_h": output_step_h, "max_step_h": max_step_h, "rtol": 1e-7, "atol": 1e-9}
    result = simulate(
        model, initial_cfu_ml, regimen, duration_days=duration_days,
        snapshot_days=snapshot_days, others_cfu_ml=others_cfu_ml, **solver,
    )
    save_results(result, model, patient, initial_cfu_ml, regimen, snapshot_days,
                 detection_limit, output_dir, solver, initial_provenance)
    return result


def run_example(
    output_dir: Path = Path("eredmenyek"), regimen: Sequence[DailyDose] | None = None,
    snapshot_days: Sequence[float] | None = None,
) -> Result:
    if parameters.INITIAL_MODE == "relative_anchor":
        initial, others = initial_from_relative(
            parameters.INITIAL_PERCENTAGES, parameters.ANCHOR_SPECIES, parameters.ANCHOR_LOG10_CFU_ML,
        )
        provenance = "Estimated CFU/mL from hypothetical absolute anchor and relative proportions."
    elif parameters.INITIAL_MODE == "log10":
        initial = initial_from_log10(parameters.INITIAL_LOG10_CFU_ML)
        others = parameters.OTHERS_CFU_ML
        provenance = "User-provided log10 CFU/mL; example values are hypothetical."
    elif parameters.INITIAL_MODE == "cfu":
        initial = parameters.INITIAL_CFU_ML
        others = parameters.OTHERS_CFU_ML
        provenance = "User-provided CFU/mL; example values are hypothetical."
    else:
        raise ValueError("INITIAL_MODE must be relative_anchor, log10, or cfu.")
    return run_simulation(
        parameters.build_model(), initial, parameters.REGIMEN if regimen is None else regimen,
        patient=parameters.PATIENT, others_cfu_ml=others, duration_days=parameters.DURATION_DAYS,
        snapshot_days=parameters.SNAPSHOT_DAYS if snapshot_days is None else snapshot_days,
        output_dir=output_dir, output_step_h=parameters.OUTPUT_STEP_H, max_step_h=parameters.MAX_STEP_H,
        detection_limit=parameters.DETECTION_LIMIT_CFU_ML, initial_provenance=provenance,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "eredmenyek")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--regimen-csv", type=Path)
    group.add_argument("--no-antibiotics", action="store_true")
    parser.add_argument("--snapshots", type=float, nargs="+")
    args = parser.parse_args()
    regimen = None
    if args.regimen_csv is not None:
        regimen = load_regimen_csv(args.regimen_csv)
    elif args.no_antibiotics:
        regimen = ()
    result = run_example(args.output, regimen, args.snapshots)
    print(f"Saved {len(result.time_h)} samples, {len(result.events)} doses to {args.output.resolve()}")
    print(f"Final relative abundance sum: {result.percentages[-1].sum():.12g}%")


if __name__ == "__main__":
    main()
