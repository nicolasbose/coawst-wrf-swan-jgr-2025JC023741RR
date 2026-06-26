"""Drag-coefficient plots using the precomputed OWF drag NetCDF.

Reads drag_owf_JFM2024.nc and generates the
figures to the OWF energy-map boundary.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Drag-coefficient plots from precomputed OWF drag NetCDF."
    )
    parser.add_argument(
        "--drag-nc",
        default="database/drag_owf_JFM2024.nc",
        help="Precomputed drag/wave-age NetCDF produced by drag_plots.py.",
    )
    parser.add_argument(
        "--output-dir",
        default="figures",
        help="Directory where the figures will be saved.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Binning helpers
# ---------------------------------------------------------------------------

def separete_variable_wind_bins(
    wspd: np.ndarray,
    main_data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    wind_bins = np.arange(1, 13, 0.5)
    bin_centers = (wind_bins[:-1] + wind_bins[1:]) / 2.0
    means: list[float] = []
    stds: list[float] = []
    for left, right in zip(wind_bins[:-1], wind_bins[1:]):
        mask = (wspd >= left) & (wspd < right)
        values = main_data[mask]
        if values.size:
            means.append(float(np.nanmean(values)))
            stds.append(float(np.nanstd(values)))
        else:
            means.append(np.nan)
            stds.append(np.nan)
    return bin_centers, np.asarray(means), np.asarray(stds)


def bin_wind_wave_alignment(angle_wind_wave: np.ndarray) -> np.ndarray:
    bins = np.array([0.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360.0])
    bin_labels = np.array([0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 150.0, 120.0, 90.0, 60.0, 30.0, 0.0])
    clipped = np.mod(angle_wind_wave, 360.0)
    bin_indices = np.digitize(clipped, bins, right=True)
    return np.array([bin_labels[index - 1] for index in bin_indices], dtype=float)


def reorder_by_alignment_bins(
    angle_raw: np.ndarray,
    binned_directions: np.ndarray,
    *arrays: np.ndarray,
) -> tuple[np.ndarray, ...]:
    ranges = [
        (0.0, 30.0, False),
        (30.0, 60.0, False),
        (60.0, 90.0, False),
        (90.0, 120.0, False),
        (120.0, 150.0, False),
        (150.0, 180.0, True),
    ]
    grouped: list[list[np.ndarray]] = [[] for _ in range(len(arrays) + 1)]
    for lower, upper, include_upper in ranges:
        if include_upper:
            mask = (binned_directions >= lower) & (binned_directions <= upper)
        else:
            mask = (binned_directions >= lower) & (binned_directions < upper)
        grouped[0].append(angle_raw[mask])
        for idx, values in enumerate(arrays, start=1):
            grouped[idx].append(values[mask])
    return tuple(np.concatenate(parts) for parts in grouped)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_drag_netcdf(path: str | Path) -> dict[str, np.ndarray | None]:
    ds = xr.open_dataset(path)

    wspd10_wrfswan  = ds.wspd10_wrfswan.values.flatten()
    wspd10_wrfstand = ds.wspd10_wrfstand.values.flatten()
    wspd100_wrfswan  = ds.wspd100_wrfswan.values.flatten()
    wspd100_wrfstand = ds.wspd100_wrfstand.values.flatten()
    cd_wrfswan  = ds.cd_wrfswan.values.flatten()
    cd_wrfstand = ds.cd_wrfstand.values.flatten()
    d_u10     = ds.d_u10.values.flatten()
    d_wspd100 = ds.d_wspd100.values.flatten()
    wave_age  = ds.wave_age.values.flatten()

    valid = (
        np.isfinite(wspd10_wrfswan) & np.isfinite(wspd10_wrfstand)
        & np.isfinite(wspd100_wrfswan) & np.isfinite(wspd100_wrfstand)
        & np.isfinite(cd_wrfswan) & np.isfinite(cd_wrfstand)
        & np.isfinite(wave_age)
    )

    angle_wind_wave = None
    angle_wind_wave_binned = None
    if "angle_wind_wave" in ds:
        angle_raw = ds.angle_wind_wave.values.flatten()
        valid = valid & np.isfinite(angle_raw)
        angle_wind_wave = angle_raw
        angle_wind_wave_binned = bin_wind_wave_alignment(angle_raw)

    data: dict[str, np.ndarray | None] = {
        "wspd10_wrfswan":  wspd10_wrfswan[valid],
        "wspd10_wrfstand": wspd10_wrfstand[valid],
        "wspd100_wrfswan":  wspd100_wrfswan[valid],
        "wspd100_wrfstand": wspd100_wrfstand[valid],
        "cd_wrfswan":  cd_wrfswan[valid],
        "cd_wrfstand": cd_wrfstand[valid],
        "d_u10":     d_u10[valid],
        "d_wspd100": d_wspd100[valid],
        "wave_age": wave_age[valid],
        "angle_wind_wave": angle_wind_wave[valid] if angle_wind_wave is not None else None,
        "angle_wind_wave_binned": (
            angle_wind_wave_binned[valid] if angle_wind_wave_binned is not None else None
        ),
    }

    (
        data["bin_centers_wrfswan"],
        data["mean_wrfswan"],
        data["std_dev_wrfswan"],
    ) = separete_variable_wind_bins(data["wspd10_wrfswan"], data["cd_wrfswan"] * 1000.0)
    (
        data["bin_centers_wrfstand"],
        data["mean_wrfstand"],
        data["std_dev_wrfstand"],
    ) = separete_variable_wind_bins(data["wspd10_wrfstand"], data["cd_wrfstand"] * 1000.0)

    return data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_main_drag_figure(
    data: dict[str, np.ndarray | None],
    output_path: str | Path,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mask = (
        (data["wspd10_wrfswan"] >= 1.0)
        & (data["wspd10_wrfswan"] <= 13.0)
        & (data["wspd10_wrfstand"] >= 1.0)
        & (data["wspd10_wrfstand"] <= 13.0)
    )
    d_u10 = data["d_u10"][mask]
    lim = float(np.nanpercentile(np.abs(d_u10), 98))
    vmax = max(0.5, lim)
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    cmap = plt.cm.coolwarm

    plt.rcParams.update({"font.size": 16})
    fig, ax = plt.subplots(1, 2, figsize=(25, 8), sharex=True, sharey=True)

    ax[0].scatter(
        data["wspd10_wrfswan"][mask],
        data["cd_wrfswan"][mask] * 1000.0,
        c=d_u10,
        cmap=cmap,
        norm=norm,
        s=18,
        linewidths=0,
    )
    ax[0].errorbar(
        data["bin_centers_wrfswan"],
        data["mean_wrfswan"],
        yerr=data["std_dev_wrfswan"],
        color="k",
        fmt="o--",
        capsize=5,
        label="Mean ± SD",
    )
    ax[0].set_title("Coupled: WRF-SWAN", pad=10)
    ax[0].set_xlabel(r"Wspd$_{10}$ [m s$^{-1}$]")
    ax[0].set_ylabel(r"1000 × C$_{D}$")
    ax[0].grid(True)
    ax[0].legend(loc="upper left")

    ax[1].scatter(
        data["wspd10_wrfstand"][mask],
        data["cd_wrfstand"][mask] * 1000.0,
        c=d_u10,
        cmap=cmap,
        norm=norm,
        s=18,
        linewidths=0,
    )
    ax[1].errorbar(
        data["bin_centers_wrfstand"],
        data["mean_wrfstand"],
        yerr=data["std_dev_wrfstand"],
        color="k",
        fmt="o--",
        capsize=5,
        label="Mean ± SD",
    )
    ax[1].set_title("Stand-alone: WRF", pad=10)
    ax[1].set_xlabel(r"Wspd$_{10}$ [m s$^{-1}$]")
    ax[1].set_ylabel(r"1000 × C$_{D}$")
    ax[1].grid(True)
    ax[1].legend(loc="upper left")

    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        location="right",
        fraction=0.05,
        pad=0.02,
    )
    cbar.set_label(
        r"$\Delta$ U$_{10}$ = Wspd$_{10}$(WRF-SWAN) − Wspd$_{10}$(WRF) [m s$^{-1}$]"
    )

    fig.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close(fig)
    return output_path


def make_wave_age_drag_figure(
    data: dict[str, np.ndarray | None],
    output_path: str | Path,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mask = (
        (data["wspd10_wrfswan"] >= 1.0)
        & (data["wspd10_wrfswan"] <= 13.0)
        & (data["wspd10_wrfstand"] >= 1.0)
        & (data["wspd10_wrfstand"] <= 13.0)
    )
    wave_age = data["wave_age"][mask]
    vmin = float(np.nanmin(wave_age))
    vmax = float(np.nanmax(wave_age))
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.025, vmax=vmax)

    if data["angle_wind_wave"] is not None and data["angle_wind_wave_binned"] is not None:
        (
            ordered_angle,
            ordered_wspd10,
            ordered_cd,
            ordered_wave_age,
        ) = reorder_by_alignment_bins(
            data["angle_wind_wave"][mask],
            data["angle_wind_wave_binned"][mask],
            data["wspd10_wrfswan"][mask],
            data["cd_wrfswan"][mask],
            wave_age,
        )
    else:
        ordered_angle = None
        ordered_wspd10 = data["wspd10_wrfswan"][mask]
        ordered_cd = data["cd_wrfswan"][mask]
        ordered_wave_age = wave_age

    plt.rcParams.update({"font.size": 16})
    if data["angle_wind_wave"] is not None:
        fig, ax = plt.subplots(1, 2, figsize=(25, 8), sharex=True, sharey=True)
        axes = np.asarray(ax)
    else:
        fig, ax = plt.subplots(figsize=(12, 7))
        axes = np.asarray([ax])

    axes[0].scatter(
        ordered_wspd10,
        ordered_cd * 1000.0,
        c=ordered_wave_age,
        cmap="RdYlBu_r",
        norm=norm,
        s=22,
        linewidths=0,
        alpha=0.9,
    )
    axes[0].errorbar(
        data["bin_centers_wrfswan"],
        data["mean_wrfswan"],
        yerr=data["std_dev_wrfswan"],
        color="k",
        fmt="o--",
        capsize=5,
        label="Mean ± Std Dev",
    )
    axes[0].set_title("Coupled: WRF-SWAN", pad=8)
    axes[0].set_xlabel(r"Wspd$_{10}$ [m s$^{-1}$]")
    axes[0].set_ylabel(r"1000 × C$_{D}$")
    axes[0].grid(True, alpha=0.6)
    axes[0].legend(loc="upper left")
    axes[0].set_xlim(0.5, 12.6)
    axes[0].set_ylim(0.68, 2.76)

    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap="RdYlBu_r"),
        ax=axes[0],
        location="right",
        fraction=0.05,
        pad=0.02,
    )
    cbar.set_label(r"$u_{*}/C_{p}$")

    if data["angle_wind_wave"] is not None:
        norm_angle = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=90.0, vmax=180.0)
        axes[1].scatter(
            ordered_wspd10,
            ordered_cd * 1000.0,
            c=ordered_angle,
            cmap="coolwarm",
            norm=norm_angle,
            s=22,
            linewidths=0,
            alpha=0.9,
        )
        axes[1].errorbar(
            data["bin_centers_wrfswan"],
            data["mean_wrfswan"],
            yerr=data["std_dev_wrfswan"],
            color="k",
            fmt="o--",
            capsize=5,
            label="Mean ± Std Dev",
        )
        axes[1].set_title("Coupled: WRF-SWAN", pad=8)
        axes[1].set_xlabel(r"Wspd$_{10}$ [m s$^{-1}$]")
        axes[1].set_ylabel(r"1000 × C$_{D}$")
        axes[1].grid(True, alpha=0.6)
        axes[1].legend(loc="upper left")
        axes[1].set_xlim(0.5, 12.6)
        axes[1].set_ylim(0.68, 2.76)

        cbar2 = fig.colorbar(
            plt.cm.ScalarMappable(norm=norm_angle, cmap="coolwarm"),
            ax=axes[1],
            location="right",
            fraction=0.05,
            pad=0.02,
        )
        cbar2.set_label("Wind-Wave alignment [°]")

    fig.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> None:
    args = parse_args()

    print(f"Loading precomputed drag NetCDF: {args.drag_nc}")
    data = load_drag_netcdf(path=args.drag_nc)

    output_dir = Path(args.output_dir)
    main_path = make_main_drag_figure(
        data=data,
        output_path=output_dir / "Figure4_CDN_vs_U10N_WRF-SWAN_vs_WRF_owf.png",
    )
    wave_age_path = make_wave_age_drag_figure(
        data=data,
        output_path=output_dir / "Cd10N_wave_age_wrfswan_owf.png",
    )

    print(f"Saved main drag figure to: {main_path.resolve()}")
    print(f"Saved wave-age drag figure to: {wave_age_path.resolve()}")
    print(f"Samples used: {data['wspd10_wrfswan'].size}")
    print(
        f"Coupled Wspd10 range: "
        f"{np.nanmin(data['wspd10_wrfswan']):.4f} to {np.nanmax(data['wspd10_wrfswan']):.4f}"
    )
    print(
        f"Stand-alone Wspd10 range: "
        f"{np.nanmin(data['wspd10_wrfstand']):.4f} to {np.nanmax(data['wspd10_wrfstand']):.4f}"
    )
    print(
        f"Coupled Cd range: "
        f"{np.nanmin(data['cd_wrfswan'] * 1000.0):.4f} to "
        f"{np.nanmax(data['cd_wrfswan'] * 1000.0):.4f}"
    )
    print(
        f"Stand-alone Cd range: "
        f"{np.nanmin(data['cd_wrfstand'] * 1000.0):.4f} to "
        f"{np.nanmax(data['cd_wrfstand'] * 1000.0):.4f}"
    )


if __name__ == "__main__":
    main()
