"""Standalone wave-stress plots for the JGR submission package.

This script reads the prepared database files,
derives the inverse wave age from the coupled WRF friction velocity and the
peak period inferred from the 1D spectrum, and reproduces the wave-stress
scatter plots with the 100 m wind-speed difference used as the color scale.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import xarray as xr

G = 9.81


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create standalone wave-stress plots from the database files."
    )
    parser.add_argument(
        "--wind-nc",
        default="database/wrf_wind_owf_JFM2024.nc",
        help="Precomputed WRF wind variables at OWF coordinates (from extract_wrf_wind_owf.py).",
    )
    parser.add_argument(
        "--spec",
        default="database/spec_1d_tau_with_tau_chen_owf.nc",
        help="SWAN post-process spectrum file with tm01 and tau_chen.",
    )
    parser.add_argument(
        "--output-dir",
        default="figures",
        help="Directory where the figure will be saved.",
    )
    parser.add_argument(
        "--wave-direction-var",
        default="wave_dir",
        help="Wave-direction variable in the SWAN post-process file for the angle figure.",
    )
    return parser.parse_args()


def estimate_wave_age(u: np.ndarray, tp: np.ndarray) -> np.ndarray:
    cp = (G * tp) / (2.0 * np.pi)
    return u / cp


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


def compute_wind_wave_alignment(
    wind_direction: np.ndarray,
    wave_direction: np.ndarray,
) -> np.ndarray:
    angle_wind_wave: list[float] = []
    for wind_dir, wave_dir in zip(wind_direction, wave_direction):
        theta1 = np.deg2rad(wind_dir)
        theta2 = np.deg2rad(wave_dir)
        u1, v1 = np.cos(theta1), np.sin(theta1)
        u2, v2 = np.cos(theta2), np.sin(theta2)
        cos_theta = np.dot(u1, u2) + np.dot(v1, v2)
        norm1 = np.sqrt(np.dot(u1, u1) + np.dot(v1, v1))
        norm2 = np.sqrt(np.dot(u2, u2) + np.dot(v2, v2))
        cos_similarity = cos_theta / (norm1 * norm2)
        cos_similarity = np.clip(cos_similarity, -1.0, 1.0)
        angle_degrees = np.arccos(cos_similarity) * 180.0 / np.pi
        angle_wind_wave.append(float(angle_degrees))
    return np.asarray(angle_wind_wave)


def build_wave_stress_arrays_notebook_exact(
    wind_nc_path: str | Path,
    spec_path: str | Path,
    wave_direction_var: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray | None,
]:
    wind = xr.open_dataset(wind_nc_path)
    spec = xr.open_dataset(spec_path)

    common_times = np.intersect1d(wind.time.values, spec.time.values)
    wind = wind.sel(time=common_times)
    spec = spec.sel(time=common_times)

    ust_wrfswan     = wind.ust_wrfswan.values.flatten()
    wspd100_wrfswan  = wind.wspd100_wrfswan.values.flatten()
    wspd100_wrfstand = wind.wspd100_wrfstand.values.flatten()
    wind_direction   = wind.wdir_wrfswan.values.flatten()

    wave_age     = estimate_wave_age(ust_wrfswan, spec.tm01.values.flatten())
    wspd100_diff = wspd100_wrfswan - wspd100_wrfstand
    tau_wave     = spec.tau_chen.values.flatten()

    if wave_direction_var in spec:
        wave_direction = spec[wave_direction_var].values.flatten()
    else:
        wave_direction = None

    valid = (
        np.isfinite(wave_age)
        & np.isfinite(tau_wave)
        & np.isfinite(wspd100_wrfswan)
        & np.isfinite(wspd100_diff)
    )
    if wave_direction is not None:
        valid = valid & np.isfinite(wind_direction) & np.isfinite(wave_direction)

    wave_age = wave_age[valid]
    tau_wave = tau_wave[valid]
    wspd100_wrfswan = wspd100_wrfswan[valid]
    wspd100_diff = wspd100_diff[valid]
    if wave_direction is not None:
        angle_wind_wave = compute_wind_wave_alignment(
            wind_direction[valid],
            wave_direction[valid],
        )
    else:
        angle_wind_wave = None

    negative_mask = wspd100_diff < 0
    bin_centers_wrfswan, mean_wrfswan, std_dev_wrfswan = separete_variable_wind_bins(
        wave_age[negative_mask],
        tau_wave[negative_mask],
    )
    positive_mask = wspd100_diff > 0
    (
        bin_centers_wrfswan_2,
        mean_wrfswan_2,
        std_dev_wrfswan_2,
    ) = separete_variable_wind_bins(
        wave_age[positive_mask],
        tau_wave[positive_mask],
    )

    return (
        wave_age,
        tau_wave,
        wspd100_wrfswan,
        wspd100_diff,
        bin_centers_wrfswan,
        mean_wrfswan,
        std_dev_wrfswan,
        bin_centers_wrfswan_2,
        mean_wrfswan_2,
        std_dev_wrfswan_2,
        angle_wind_wave,
    )


def make_wave_stress_figure_notebook_exact(
    wave_age: np.ndarray,
    tau_wave: np.ndarray,
    wspd100: np.ndarray,
    delta_wspd100: np.ndarray,
    bin_centers_wrfswan: np.ndarray,
    mean_wrfswan: np.ndarray,
    std_dev_wrfswan: np.ndarray,
    bin_centers_wrfswan_2: np.ndarray,
    mean_wrfswan_2: np.ndarray,
    std_dev_wrfswan_2: np.ndarray,
    output_path: str | Path,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.size": 14})
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(15, 8),
        sharex=False,
        sharey=True,
    )

    scatter_top = axes[0].scatter(
        wave_age,
        tau_wave,
        c=wspd100,
        cmap="RdYlBu_r",
        s=22,
    )
    axes[0].axhline(0.0, color="black", linestyle="--", linewidth=1)
    axes[0].set_ylabel(r"$\tau_{wave(0)} / u_{*}^{2}$")
    axes[0].grid(True, alpha=0.3)

    norm = mcolors.TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
    negative_mask = delta_wspd100 < 0
    axes[1].errorbar(
        bin_centers_wrfswan,
        mean_wrfswan,
        yerr=std_dev_wrfswan,
        color="k",
        fmt="o--",
        capsize=5,
        label="Mean ± Std Dev",
    )
    axes[1].errorbar(
        bin_centers_wrfswan_2,
        mean_wrfswan_2,
        yerr=std_dev_wrfswan_2,
        color="k",
        fmt="o--",
        capsize=5,
        label="Mean ± Std Dev",
    )
    scatter_bottom = axes[1].scatter(
        wave_age,
        tau_wave,
        c=delta_wspd100,
        cmap="coolwarm",
        norm=norm,
        s=22,
    )
    axes[1].axhline(0.0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel(r"$u_{*} / C_{p}$")
    axes[1].set_ylabel(r"$\tau_{wave(0)} / u_{*}^{2}$")
    axes[1].grid(True, alpha=0.3)
    cbar_top = fig.colorbar(scatter_top, ax=axes[0], location="right", fraction=0.15, pad=0.01)
    cbar_top.set_label(r"Wspd$_{100}$ [m.$s^{-1}$]")

    cbar_bottom = fig.colorbar(
        scatter_bottom,
        ax=axes[1],
        location="right",
        fraction=0.15,
        pad=0.01,
    )
    cbar_bottom.set_label(r"$\Delta$ Wspd$_{100}$ [m.$s^{-1}]$")

    fig.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close(fig)
    return output_path


def make_wind_wave_tau_angle_figure(
    wave_age: np.ndarray,
    tau_wave: np.ndarray,
    angle_wind_wave: np.ndarray,
    output_path: str | Path,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=[12, 6])
    norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=90, vmax=180)
    scatter = ax.scatter(
        wave_age,
        tau_wave,
        c=angle_wind_wave,
        cmap="coolwarm",
        norm=norm,
        s=10,
        alpha=0.9,
    )
    ax.set_ylabel(r"$\tau_{wave0} / u_{*}^{2}$")
    ax.set_xlabel(r"$u_{*} / C_{p}$")
    ax.axhline(0.0, color="black", linewidth=1, linestyle="--")

    inset_ax = inset_axes(
        ax,
        width="28%",
        height="28%",
        loc="lower right",
        bbox_to_anchor=(0.04, 0.08, 0.95, 1.15),
        bbox_transform=ax.transAxes,
    )
    inset_ax.hist(
        angle_wind_wave,
        bins=30,
        color="0.7",
        edgecolor="0.15",
        linewidth=1.0,
    )
    inset_ax.set_xlabel("Wind-Wave Alignment [$^{o}$]", fontsize=8)
    inset_ax.set_ylabel("Freq", fontsize=8)
    inset_ax.tick_params(axis="both", labelsize=7)

    fig.tight_layout()
    cbar = fig.colorbar(scatter, ax=ax, location="right", fraction=0.05, pad=0.05)
    cbar.set_label("Wind-Wave Alignment [°]")
    fig.savefig(output_path, dpi=500, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> None:
    args = parse_args()
    (
        wave_age,
        tau_wave,
        wspd100,
        delta_wspd100,
        bin_centers_wrfswan,
        mean_wrfswan,
        std_dev_wrfswan,
        bin_centers_wrfswan_2,
        mean_wrfswan_2,
        std_dev_wrfswan_2,
        angle_wind_wave,
    ) = build_wave_stress_arrays_notebook_exact(
        wind_nc_path=args.wind_nc,
        spec_path=args.spec,
        wave_direction_var=args.wave_direction_var,
    )

    output_path = Path(args.output_dir) / "tau_wind_wave_dif_100m_notebook_exact.png"
    figure_path = make_wave_stress_figure_notebook_exact(
        wave_age=wave_age,
        tau_wave=tau_wave,
        wspd100=wspd100,
        delta_wspd100=delta_wspd100,
        bin_centers_wrfswan=bin_centers_wrfswan,
        mean_wrfswan=mean_wrfswan,
        std_dev_wrfswan=std_dev_wrfswan,
        bin_centers_wrfswan_2=bin_centers_wrfswan_2,
        mean_wrfswan_2=mean_wrfswan_2,
        std_dev_wrfswan_2=std_dev_wrfswan_2,
        output_path=output_path,
    )
    if angle_wind_wave is not None:
        angle_output_path = Path(args.output_dir) / "wind_wave_tau_angle_100m.png"
        angle_figure_path = make_wind_wave_tau_angle_figure(
            wave_age=wave_age,
            tau_wave=tau_wave,
            angle_wind_wave=angle_wind_wave,
            output_path=angle_output_path,
        )
        print(f"Saved angle figure to: {angle_figure_path.resolve()}")
    else:
        print(
            "Angle figure skipped: "
            f"wave-direction variable '{args.wave_direction_var}' was not found in the spec file."
        )

    print(f"Saved figure to: {figure_path.resolve()}")
    print(f"Samples used: {wave_age.size}")
    print(f"Wind file: {Path(args.wind_nc).resolve()}")
    print(f"Spec file: {Path(args.spec).resolve()}")
    print(
        "Inverse wave age range: "
        f"{np.nanmin(wave_age):.4f} to {np.nanmax(wave_age):.4f}"
    )
    print(f"Tau-wave range: {np.nanmin(tau_wave):.4f} to {np.nanmax(tau_wave):.4f}")
    print(f"Wspd100 range: {np.nanmin(wspd100):.4f} to {np.nanmax(wspd100):.4f}")
    print(
        "Delta Wspd100 range: "
        f"{np.nanmin(delta_wspd100):.4f} to {np.nanmax(delta_wspd100):.4f}"
    )


if __name__ == "__main__":
    main()
