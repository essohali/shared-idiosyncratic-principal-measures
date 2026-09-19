#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FINAL extraction script for the SJS application
================================================

Allen Visual Coding Neuropixels / DANDI
Dandiset : 000022
Version  : 0.251116.2247
Session  : 794812542
Asset    : sub-774672354/sub-774672354_ses-794812542.nwb

This script DOES NOT download the complete NWB file. It streams the remote
asset and exports only the data required for the paper.

Final application design
------------------------
Five areas actually available in this NWB:
    VISp, VISl, VISrl, VISal, VISam

Quality rule:
    amplitude_cutoff < 0.1
    presence_ratio > 0.9
    isi_violations < 0.5

Balanced population:
    40 quality units per area, selected deterministically by unit_id.
    This is independent of stimulus response/covariance and avoids cherry-picking.

Stimulus:
    interval table = drifting_gratings_75_repeats_presentations
    orientation = 0 degrees
    contrast = 0.8
    temporal frequency = 2 Hz
    spatial frequency = 0.04 cycles/degree
    expected 75 repeated presentations
    analysis window = [0, 2] seconds relative to onset

Outputs:
    Allen_794812542_application.xlsx
    units_quality.csv
    selected_units.csv
    trials_75.csv
    spikes_application.csv
    metadata.csv
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# FIXED ANALYSIS SETTINGS
# ---------------------------------------------------------------------
DANDISET_ID = "000022"
DANDISET_VERSION = "0.251116.2247"
SESSION_ID = 794812542
ASSET_PATH = "sub-774672354/sub-774672354_ses-794812542.nwb"

AREAS = ["VISp", "VISl", "VISrl", "VISal", "VISam"]
N_UNITS_PER_AREA = 40

STIM_TABLE = "drifting_gratings_75_repeats_presentations"
ORIENTATION_DEG = 0.0
WINDOW_START = 0.0
WINDOW_END = 2.0

OUTDIR = Path("Allen_794812542_final_data")
OUTDIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# DEPENDENCIES
# ---------------------------------------------------------------------
def load_dependencies():
    try:
        from dandi.dandiapi import DandiAPIClient
        import remfile
        import h5py
        from pynwb import NWBHDF5IO
    except ImportError as exc:
        raise SystemExit(
            "\nMissing dependency.\n"
            "Run this in Jupyter:\n\n"
            "%pip install dandi pynwb h5py remfile pandas numpy openpyxl\n\n"
            "Then restart the kernel and rerun this script.\n"
        ) from exc

    return DandiAPIClient, remfile, h5py, NWBHDF5IO


# ---------------------------------------------------------------------
# REMOTE NWB
# ---------------------------------------------------------------------
def locate_asset(DandiAPIClient):
    print("1) Locating published DANDI asset...")
    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(DANDISET_ID, DANDISET_VERSION)
        asset = dandiset.get_asset_by_path(ASSET_PATH)

        url = asset.get_content_url(
            follow_redirects=1,
            strip_query=True
        )

        asset_id = (
            getattr(asset, "identifier", None)
            or getattr(asset, "asset_id", None)
        )
        size = getattr(asset, "size", None)

    print(f"   Dandiset: {DANDISET_ID}, version {DANDISET_VERSION}")
    print(f"   Session: {SESSION_ID}")
    if size:
        print(f"   Remote NWB size: {size / 1024**3:.2f} GiB")
    return url, asset_id, size


def open_remote_nwb(url, remfile, h5py, NWBHDF5IO):
    print("2) Opening NWB remotely; no full NWB download...")
    rf = remfile.File(url)
    hf = h5py.File(rf, "r")
    io = NWBHDF5IO(file=hf, load_namespaces=True)
    nwb = io.read()
    return rf, hf, io, nwb


# ---------------------------------------------------------------------
# UNITS / ANATOMY
# ---------------------------------------------------------------------
def normalize_area(value):
    if pd.isna(value):
        return None
    s = str(value).strip()
    for area in AREAS:
        if s == area or area in s:
            return area
    return s


def build_unit_metadata(units_table):
    ids = np.asarray(units_table.id[:])
    d = {"unit_id": ids}

    # Exclude spike_times here: it is a large ragged array read later only
    # for the 200 selected units.
    for col in units_table.colnames:
        if col == "spike_times":
            continue
        try:
            values = units_table[col][:]
            if len(values) == len(ids):
                # Avoid bringing large ragged waveform/spike amplitude arrays
                # into the compact metadata if they are non-scalar per unit.
                scalar_like = True
                sample = values[: min(5, len(values))]
                for x in sample:
                    if isinstance(x, (list, tuple, np.ndarray)) and np.ndim(x) > 0:
                        scalar_like = False
                        break
                if scalar_like:
                    d[col] = list(values)
        except Exception:
            pass

    return pd.DataFrame(d)


def read_units(nwb):
    print("3) Reading units and anatomical locations...")

    if nwb.units is None:
        raise RuntimeError("No units table found.")

    units_table = nwb.units
    units = build_unit_metadata(units_table)

    if "peak_channel_id" not in units.columns:
        raise RuntimeError(
            "peak_channel_id is absent from the units table."
        )

    if nwb.electrodes is None:
        raise RuntimeError("No electrodes table found.")

    electrodes = nwb.electrodes.to_dataframe().reset_index()

    id_col = next(
        (c for c in ["id", "channel_id", "ecephys_channel_id", "index"]
         if c in electrodes.columns),
        None
    )
    structure_col = next(
        (c for c in ["ecephys_structure_acronym",
                     "structure_acronym",
                     "location"]
         if c in electrodes.columns),
        None
    )

    if id_col is None or structure_col is None:
        raise RuntimeError(
            "Cannot link units to anatomy.\n"
            f"Electrode columns: {list(electrodes.columns)}"
        )

    anatomy = (
        electrodes[[id_col, structure_col]]
        .drop_duplicates()
        .rename(columns={
            id_col: "peak_channel_id",
            structure_col: "area"
        })
    )

    units = units.merge(
        anatomy,
        on="peak_channel_id",
        how="left"
    )
    units["area"] = units["area"].map(normalize_area)
    units = units[units["area"].isin(AREAS)].copy()

    # Strict prespecified quality rule.
    required = ["amplitude_cutoff", "presence_ratio", "isi_violations"]
    missing = [c for c in required if c not in units.columns]
    if missing:
        raise RuntimeError(
            "Cannot apply the prespecified quality rule because columns "
            f"are missing: {missing}"
        )

    amp = pd.to_numeric(units["amplitude_cutoff"], errors="coerce")
    presence = pd.to_numeric(units["presence_ratio"], errors="coerce")
    isi = pd.to_numeric(units["isi_violations"], errors="coerce")

    units["passes_quality"] = (
        amp.notna()
        & presence.notna()
        & isi.notna()
        & (amp < 0.1)
        & (presence > 0.9)
        & (isi < 0.5)
    )

    units_quality = (
        units.loc[units["passes_quality"]]
        .copy()
        .sort_values(["area", "unit_id"])
        .reset_index(drop=True)
    )

    counts = (
        units.groupby("area")["unit_id"].size()
        .reindex(AREAS, fill_value=0)
        .rename("n_units")
        .to_frame()
    )
    counts["n_quality_units"] = (
        units_quality.groupby("area")["unit_id"].size()
        .reindex(AREAS, fill_value=0)
    )
    counts = counts.reset_index()

    print(counts.to_string(index=False))

    insufficient = counts.loc[
        counts["n_quality_units"] < N_UNITS_PER_AREA
    ]
    if len(insufficient):
        raise RuntimeError(
            "\nAt least one area has fewer than "
            f"{N_UNITS_PER_AREA} quality units:\n"
            + insufficient.to_string(index=False)
        )

    # Deterministic balanced selection, independent of response.
    selected_units = (
        units_quality
        .sort_values(["area", "unit_id"])
        .groupby("area", group_keys=False)
        .head(N_UNITS_PER_AREA)
        .copy()
        .reset_index(drop=True)
    )
    selected_units["selection_rule"] = (
        f"first {N_UNITS_PER_AREA} quality-passing units by ascending unit_id "
        "within area"
    )

    print(
        f"   Balanced selection: {len(selected_units)} units "
        f"({N_UNITS_PER_AREA} x {len(AREAS)} areas)."
    )

    return units_table, units_quality, selected_units, counts


# ---------------------------------------------------------------------
# STIMULUS TRIALS
# ---------------------------------------------------------------------
def table_to_df(table):
    df = table.to_dataframe().reset_index()

    if "id" in df.columns:
        df = df.rename(columns={"id": "stimulus_presentation_id"})
    elif "index" in df.columns:
        df = df.rename(columns={"index": "stimulus_presentation_id"})
    elif "stimulus_presentation_id" not in df.columns:
        df.insert(0, "stimulus_presentation_id", np.arange(len(df)))

    if "start_time" not in df.columns and "start" in df.columns:
        df = df.rename(columns={"start": "start_time"})
    if "stop_time" not in df.columns and "stop" in df.columns:
        df = df.rename(columns={"stop": "stop_time"})

    return df


def orientation_column(df):
    for col in ["orientation", "direction", "Ori"]:
        if col in df.columns:
            return col
    raise RuntimeError(
        "No orientation column found in stimulus table.\n"
        f"Columns: {list(df.columns)}"
    )


def read_trials(nwb):
    print("4) Reading 75-repeat drifting-grating presentations...")

    if nwb.intervals is None:
        raise RuntimeError("No interval tables found.")

    keys = list(nwb.intervals.keys())
    if STIM_TABLE not in keys:
        raise RuntimeError(
            f"Expected table '{STIM_TABLE}' not found.\n"
            f"Available tables: {keys}"
        )

    stim = table_to_df(nwb.intervals[STIM_TABLE])
    print("   Stimulus columns:", list(stim.columns))

    ori_col = orientation_column(stim)

    # Exact fixed stimulus condition established from the NWB table:
    # 75 repeats at orientation 0 deg and contrast 0.8.
    required = [ori_col, "contrast", "temporal_frequency", "spatial_frequency"]
    missing = [c for c in required if c not in stim.columns]
    if missing:
        raise RuntimeError(
            f"Missing columns required for the fixed stimulus condition: {missing}"
        )

    ori = pd.to_numeric(stim[ori_col], errors="coerce")
    contrast = pd.to_numeric(stim["contrast"], errors="coerce")
    tf = pd.to_numeric(stim["temporal_frequency"], errors="coerce")
    sf = pd.to_numeric(stim["spatial_frequency"], errors="coerce")

    trials = stim.loc[
        ori.notna()
        & contrast.notna()
        & tf.notna()
        & sf.notna()
        & np.isclose(ori, 0.0, atol=1e-8, rtol=0)
        & np.isclose(contrast, 0.8, atol=1e-8, rtol=0)
        & np.isclose(tf, 2.0, atol=1e-8, rtol=0)
        & np.isclose(sf, 0.04, atol=1e-8, rtol=0)
    ].copy()

    if "start_time" not in trials.columns:
        raise RuntimeError("start_time missing from stimulus table.")

    if "stop_time" in trials.columns:
        trials["stimulus_duration"] = (
            pd.to_numeric(trials["stop_time"], errors="coerce")
            - pd.to_numeric(trials["start_time"], errors="coerce")
        )

    trials = (
        trials.sort_values("start_time")
        .reset_index(drop=True)
    )
    trials["replicate_id"] = np.arange(1, len(trials) + 1)
    trials["analysis_window_start_relative"] = WINDOW_START
    trials["analysis_window_end_relative"] = WINDOW_END

    print("   Fixed stimulus condition:")
    print("      orientation        = 0 degrees")
    print("      contrast           = 0.8")
    print("      temporal frequency = 2.0 Hz")
    print("      spatial frequency  = 0.04 cycles/degree")
    print(f"   Number of selected trials: {len(trials)}")

    if len(trials) != 75:
        raise RuntimeError(
            "\nSCIENTIFIC SAFEGUARD: expected exactly 75 presentations for "
            "orientation=0 deg, contrast=0.8, temporal_frequency=2.0 Hz, "
            f"spatial_frequency=0.04 cycles/degree, but found {len(trials)}.\n"
            "The script stops instead of silently changing the fixed trial definition."
        )

    return trials, stim


# ---------------------------------------------------------------------
# SPIKES
# ---------------------------------------------------------------------
def extract_spikes(units_table, selected_units, trials):
    print("5) Streaming spike times for the 200 selected units...")

    if "spike_times" not in units_table.colnames:
        raise RuntimeError("No spike_times column in units table.")

    all_unit_ids = np.asarray(units_table.id[:]).astype(int)
    id_to_position = {
        int(uid): i for i, uid in enumerate(all_unit_ids)
    }

    unit_area = {
        int(row.unit_id): row.area
        for row in selected_units[["unit_id", "area"]]
        .itertuples(index=False)
    }

    trial_rows = list(
        trials[
            ["stimulus_presentation_id", "replicate_id", "start_time"]
        ].itertuples(index=False, name=None)
    )

    rows = []
    n_units = len(selected_units)

    for j, uid in enumerate(
        selected_units["unit_id"].astype(int),
        start=1
    ):
        pos = id_to_position.get(int(uid))
        if pos is None:
            raise RuntimeError(f"Selected unit {uid} not found in NWB units table.")

        spike_times = np.asarray(
            units_table["spike_times"][pos],
            dtype=float
        )

        for presentation_id, replicate_id, onset in trial_rows:
            onset = float(onset)
            left = np.searchsorted(
                spike_times,
                onset + WINDOW_START,
                side="left"
            )
            right = np.searchsorted(
                spike_times,
                onset + WINDOW_END,
                side="right"
            )

            ss = spike_times[left:right]
            if len(ss):
                relative = ss - onset
                rows.extend([
                    (
                        int(replicate_id),
                        presentation_id,
                        int(uid),
                        unit_area[int(uid)],
                        float(abs_t),
                        float(rel_t),
                    )
                    for abs_t, rel_t in zip(ss, relative)
                ])

        if j % 20 == 0 or j == n_units:
            print(f"   Units processed: {j}/{n_units}")

    spikes = pd.DataFrame(
        rows,
        columns=[
            "replicate_id",
            "stimulus_presentation_id",
            "unit_id",
            "area",
            "spike_time_absolute",
            "spike_time_relative",
        ],
    )

    if len(spikes):
        spikes = (
            spikes.sort_values([
                "replicate_id",
                "area",
                "unit_id",
                "spike_time_relative",
            ])
            .reset_index(drop=True)
        )

    print(f"   Extracted spike-event rows: {len(spikes):,}")
    return spikes


# ---------------------------------------------------------------------
# METADATA / EXPORT
# ---------------------------------------------------------------------
def unique_summary(df, col):
    if col not in df.columns:
        return ""
    vals = df[col].dropna().unique().tolist()
    return "; ".join(map(str, vals[:20]))


def build_metadata(asset_id, asset_size, counts, trials, spikes):
    metadata = [
        ("dataset", "Allen Brain Observatory Visual Coding Neuropixels"),
        ("dandiset_id", DANDISET_ID),
        ("dandiset_version", DANDISET_VERSION),
        ("session_id", SESSION_ID),
        ("asset_path", ASSET_PATH),
        ("asset_id", asset_id or ""),
        ("remote_asset_size_bytes", asset_size or ""),
        ("access_method", "remote DANDI/S3 NWB streaming; complete NWB not downloaded"),
        ("areas", ", ".join(AREAS)),
        ("q", len(AREAS)),
        ("quality_rule",
         "amplitude_cutoff < 0.1; presence_ratio > 0.9; isi_violations < 0.5"),
        ("unit_selection_rule",
         f"{N_UNITS_PER_AREA} quality units per area; ascending unit_id; "
         "selection independent of stimulus response/covariance"),
        ("n_units_per_area", N_UNITS_PER_AREA),
        ("n_selected_units_total", N_UNITS_PER_AREA * len(AREAS)),
        ("stimulus_interval_table", STIM_TABLE),
        ("orientation_deg", ORIENTATION_DEG),
        ("n_trials", len(trials)),
        ("analysis_window_seconds", f"[{WINDOW_START}, {WINDOW_END}]"),
        ("n_spike_event_rows", len(spikes)),
        ("contrast_values_in_selected_trials", unique_summary(trials, "contrast")),
        ("temporal_frequency_values_in_selected_trials",
         unique_summary(trials, "temporal_frequency")),
        ("spatial_frequency_values_in_selected_trials",
         unique_summary(trials, "spatial_frequency")),
        ("scientific_note",
         "The five-area design reflects actual anatomical coverage in this NWB; "
         "VISpm was absent and is not included."),
    ]

    # Preserve observed quality counts in metadata.
    for row in counts.itertuples(index=False):
        metadata.append(
            (f"{row.area}_quality_units_available", int(row.n_quality_units))
        )

    return pd.DataFrame(metadata, columns=["field", "value"])


def autosize_worksheet(ws, max_width=42):
    for column_cells in ws.columns:
        length = 0
        for cell in column_cells[:200]:
            try:
                length = max(length, len(str(cell.value)) if cell.value is not None else 0)
            except Exception:
                pass
        ws.column_dimensions[column_cells[0].column_letter].width = min(
            max(length + 2, 10),
            max_width
        )
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions


def export_files(
    units_quality,
    selected_units,
    trials,
    spikes,
    metadata,
):
    print("6) Writing final article dataset...")

    paths = {
        "units_quality": OUTDIR / "units_quality.csv",
        "selected_units": OUTDIR / "selected_units.csv",
        "trials": OUTDIR / "trials_75.csv",
        "spikes": OUTDIR / "spikes_application.csv",
        "metadata": OUTDIR / "metadata.csv",
        "xlsx": OUTDIR / "Allen_794812542_application.xlsx",
    }

    units_quality.to_csv(paths["units_quality"], index=False)
    selected_units.to_csv(paths["selected_units"], index=False)
    trials.to_csv(paths["trials"], index=False)
    spikes.to_csv(paths["spikes"], index=False)
    metadata.to_csv(paths["metadata"], index=False)

    # Excel is for human inspection; CSV files remain canonical for analysis.
    with pd.ExcelWriter(paths["xlsx"], engine="openpyxl") as writer:
        metadata.to_excel(writer, sheet_name="metadata", index=False)
        units_quality.to_excel(writer, sheet_name="units_quality", index=False)
        selected_units.to_excel(writer, sheet_name="selected_units", index=False)
        trials.to_excel(writer, sheet_name="trials_75", index=False)

        # Excel's hard row limit is 1,048,576. Split spike events if necessary.
        max_rows = 1_000_000
        if len(spikes) <= max_rows:
            spikes.to_excel(writer, sheet_name="spikes_application", index=False)
        else:
            for i, start in enumerate(range(0, len(spikes), max_rows), start=1):
                spikes.iloc[start:start + max_rows].to_excel(
                    writer,
                    sheet_name=f"spikes_{i:02d}",
                    index=False
                )

        # Compact data dictionary.
        dictionary = pd.DataFrame([
            ["metadata", "Provenance and frozen analysis choices."],
            ["units_quality", "All quality-passing units in the five retained visual areas."],
            ["selected_units", "Balanced deterministic sample of 40 quality units per area."],
            ["trials_75", "The 75 orientation-0 drifting-grating repeated presentations."],
            ["spikes_application", "One row per spike in [0,2] s after onset for selected units/trials."],
        ], columns=["table", "description"])
        dictionary.to_excel(writer, sheet_name="data_dictionary", index=False)

        wb = writer.book
        for ws in wb.worksheets:
            autosize_worksheet(ws)

    print("\nFinal files:")
    for p in paths.values():
        print("  ", p)

    return paths


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------
def main():
    DandiAPIClient, remfile, h5py, NWBHDF5IO = load_dependencies()

    url, asset_id, asset_size = locate_asset(DandiAPIClient)
    rf, hf, io, nwb = open_remote_nwb(
        url, remfile, h5py, NWBHDF5IO
    )

    try:
        units_table, units_quality, selected_units, counts = read_units(nwb)
        trials, _ = read_trials(nwb)
        spikes = extract_spikes(
            units_table,
            selected_units,
            trials
        )

        metadata = build_metadata(
            asset_id,
            asset_size,
            counts,
            trials,
            spikes
        )

        paths = export_files(
            units_quality,
            selected_units,
            trials,
            spikes,
            metadata,
        )

        print("\n" + "=" * 68)
        print("SUCCESS")
        print("=" * 68)
        print(f"Trials: {len(trials)}")
        print(f"Areas: {', '.join(AREAS)}")
        print(f"Selected units: {len(selected_units)}")
        print(f"Spike rows: {len(spikes):,}")
        print("\nExcel file:")
        print(paths["xlsx"])
        print(
            "\nNext: upload Allen_794812542_application.xlsx here. "
            "I can then verify the dataset and compute the article diagnostics."
        )

    finally:
        try:
            io.close()
        except Exception:
            pass
        try:
            hf.close()
        except Exception:
            pass
        try:
            rf.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
