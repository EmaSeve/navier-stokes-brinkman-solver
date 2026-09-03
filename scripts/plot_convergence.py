#!/usr/bin/env python3

"""Generate restrained spatial and temporal convergence figures."""

import argparse
import csv
import math
from html import escape
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "build" / "convergence" / "results.csv"
DEFAULT_OUTPUT_DIR = ROOT / "docs" / "convergence"

WIDTH = 1160
HEIGHT = 520
PANEL_WIDTH = 500
PLOT_WIDTH = 410
PLOT_HEIGHT = 300
PLOT_Y = 130
PANEL_X = (55, 605)

COLORS = {
    "L2_ux": "#1f77b4",
    "L2_uy": "#d55e00",
    "L2_uz": "#009e73",
    "L2_p": "#6f42c1",
}
MARKERS = {
    "L2_ux": "circle",
    "L2_uy": "square",
    "L2_uz": "diamond",
    "L2_p": "circle",
}


def load_results(results_path=RESULTS):
    try:
        with results_path.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
    except FileNotFoundError:
        raise SystemExit(
            f"results file not found: {results_path}\n"
            "Run ./scripts/run_convergence.sh first."
        ) from None
    if not rows:
        raise SystemExit(f"no convergence data in {results_path}")

    required = {"study", "N", "dt", "L2_ux", "L2_uy", "L2_uz", "L2_p"}
    missing = required.difference(rows[0])
    if missing:
        raise SystemExit(
            f"results file is missing columns: {', '.join(sorted(missing))}"
        )
    return rows


def number(row, key):
    return float(row[key])


def text(x, y, value, css_class, anchor="middle", transform=""):
    attributes = (
        f'x="{x:.1f}" y="{y:.1f}" class="{css_class}" '
        f'text-anchor="{anchor}"'
    )
    if transform:
        attributes += f' transform="{transform}"'
    return f"<text {attributes}>{escape(value)}</text>"


def superscript(exponent):
    return str(exponent).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def power_of_ten(exponent):
    return f"10{superscript(exponent)}"


def format_scale(value):
    return f"{value:.3g}"


def dash_attribute(key):
    return ' stroke-dasharray="6 3"' if key == "L2_uy" else ""


def y_exponents(low, high):
    first = math.ceil(low)
    last = math.floor(high)
    values = list(range(first, last + 1))
    if len(values) <= 6:
        return values
    stride = math.ceil(len(values) / 6)
    return values[::stride]


def marker(x, y, kind, color):
    common = f'fill="white" stroke="{color}" stroke-width="2.1"'
    if kind == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.2" {common}/>'
    if kind == "square":
        return (
            f'<rect x="{x - 4.1:.1f}" y="{y - 4.1:.1f}" width="8.2" '
            f'height="8.2" {common}/>'
        )
    points = (
        f"{x:.1f},{y - 5.0:.1f} {x + 5.0:.1f},{y:.1f} "
        f"{x:.1f},{y + 5.0:.1f} {x - 5.0:.1f},{y:.1f}"
    )
    return f'<polygon points="{points}" {common}/>'


def draw_panel(rows, series, study, panel_x, clip_id):
    scale_key = "N" if study == "spatial" else "dt"
    rows = sorted(rows, key=lambda row: number(row, scale_key))
    x_values = [number(row, scale_key) for row in rows]
    values = {
        key: [number(row, key) for row in rows]
        for key, _ in series
    }

    reference_power = -2 if study == "spatial" else 2
    reference_x = min(x_values) if study == "spatial" else max(x_values)
    coarse_index = x_values.index(reference_x)
    reference_anchor = max(data[coarse_index] for data in values.values()) * 0.55
    reference = [
        reference_anchor * (x_value / reference_x) ** reference_power
        for x_value in x_values
    ]

    all_y = [entry for data in values.values() for entry in data] + reference
    log_x_min = math.log10(min(x_values))
    log_x_max = math.log10(max(x_values))
    log_y_min = math.log10(min(all_y))
    log_y_max = math.log10(max(all_y))
    x_pad = 0.07 * (log_x_max - log_x_min)
    y_pad = 0.10 * (log_y_max - log_y_min)
    log_x_min -= x_pad
    log_x_max += x_pad
    log_y_min -= y_pad
    log_y_max += y_pad

    plot_x = panel_x + 70

    def map_x(value):
        fraction = (
            (math.log10(value) - log_x_min) / (log_x_max - log_x_min)
        )
        if study == "temporal":
            fraction = 1.0 - fraction
        return plot_x + fraction * PLOT_WIDTH

    def map_y(value):
        return PLOT_Y + PLOT_HEIGHT - (
            (math.log10(value) - log_y_min) / (log_y_max - log_y_min)
        ) * PLOT_HEIGHT

    if study == "spatial":
        title = "Spatial refinement"
        subtitle = "Δt = 10⁻⁴;  T = 10⁻³"
        x_label = "grid size  N"
        reference_label = "O(N⁻²)"
    else:
        title = "Temporal refinement"
        subtitle = "N = 256;  T = 1"
        x_label = "time step  Δt"
        reference_label = "O(Δt²)"

    output = [
        text(panel_x + PANEL_WIDTH / 2, 55, title, "panel-title"),
        text(panel_x + PANEL_WIDTH / 2, 78, subtitle, "subtitle"),
    ]

    legend_width = 105 if len(series) > 1 else 150
    legend_start = panel_x + PANEL_WIDTH / 2 - (
        len(series) * legend_width + 105
    ) / 2
    legend_y = 108
    for index, (key, label) in enumerate(series):
        x = legend_start + index * legend_width
        output.extend([
            f'<line x1="{x:.1f}" y1="{legend_y}" x2="{x + 24:.1f}" '
            f'y2="{legend_y}" stroke="{COLORS[key]}" class="series"'
            f'{dash_attribute(key)}/>',
            marker(x + 12, legend_y, MARKERS[key], COLORS[key]),
            text(x + 31, legend_y + 4, label, "legend", anchor="start"),
        ])
    ref_x = legend_start + len(series) * legend_width
    output.extend([
        f'<line x1="{ref_x:.1f}" y1="{legend_y}" x2="{ref_x + 24:.1f}" '
        f'y2="{legend_y}" class="reference"/>',
        text(ref_x + 31, legend_y + 4, reference_label, "legend", anchor="start"),
    ])

    for exponent in y_exponents(log_y_min, log_y_max):
        value = 10.0**exponent
        y = map_y(value)
        output.extend([
            f'<line x1="{plot_x}" y1="{y:.1f}" '
            f'x2="{plot_x + PLOT_WIDTH}" y2="{y:.1f}" class="grid"/>',
            text(plot_x - 10, y + 4, power_of_ten(exponent), "tick", anchor="end"),
        ])

    for value in x_values:
        x = map_x(value)
        output.extend([
            f'<line x1="{x:.1f}" y1="{PLOT_Y}" x2="{x:.1f}" '
            f'y2="{PLOT_Y + PLOT_HEIGHT}" class="grid vertical"/>',
            text(x, PLOT_Y + PLOT_HEIGHT + 23, format_scale(value), "tick"),
        ])

    output.extend([
        f'<line x1="{plot_x}" y1="{PLOT_Y + PLOT_HEIGHT}" '
        f'x2="{plot_x + PLOT_WIDTH}" y2="{PLOT_Y + PLOT_HEIGHT}" class="axis"/>',
        f'<line x1="{plot_x}" y1="{PLOT_Y}" x2="{plot_x}" '
        f'y2="{PLOT_Y + PLOT_HEIGHT}" class="axis"/>',
        f'<g clip-path="url(#{clip_id})">',
    ])

    for key, _ in series:
        points = " ".join(
            f"{map_x(x):.1f},{map_y(y):.1f}"
            for x, y in zip(x_values, values[key])
        )
        output.append(
            f'<polyline points="{points}" fill="none" stroke="{COLORS[key]}" '
            f'class="series"{dash_attribute(key)}/>'
        )
        output.extend(
            marker(map_x(x), map_y(y), MARKERS[key], COLORS[key])
            for x, y in zip(x_values, values[key])
        )

    reference_points = " ".join(
        f"{map_x(x):.1f},{map_y(y):.1f}"
        for x, y in zip(x_values, reference)
    )
    output.extend([
        f'<polyline points="{reference_points}" fill="none" class="reference"/>',
        "</g>",
        text(plot_x + PLOT_WIDTH / 2, 486, x_label, "axis-label"),
        text(
            panel_x + 17,
            PLOT_Y + PLOT_HEIGHT / 2,
            "L² error",
            "axis-label",
            transform=(
                f"rotate(-90 {panel_x + 17:.1f} "
                f"{PLOT_Y + PLOT_HEIGHT / 2:.1f})"
            ),
        ),
    ])
    return output


def write_figure(rows, output_path, series):
    spatial = [row for row in rows if row["study"] == "spatial"]
    temporal = [row for row in rows if row["study"] == "temporal"]
    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" '
        f'height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        "<defs>",
        f'<clipPath id="spatial"><rect x="{PANEL_X[0] + 70}" y="{PLOT_Y}" '
        f'width="{PLOT_WIDTH}" height="{PLOT_HEIGHT}"/></clipPath>',
        f'<clipPath id="temporal"><rect x="{PANEL_X[1] + 70}" y="{PLOT_Y}" '
        f'width="{PLOT_WIDTH}" height="{PLOT_HEIGHT}"/></clipPath>',
        "<style>",
        "text { font-family: Helvetica, Arial, sans-serif; fill: #222; }",
        ".panel-title { font-size: 17px; font-weight: 600; }",
        ".subtitle { font-size: 11.5px; fill: #555; }",
        ".legend { font-size: 11.5px; }",
        ".tick { font-size: 11px; fill: #444; }",
        ".axis-label { font-size: 13px; }",
        ".grid { stroke: #d8d8d8; stroke-width: 1; }",
        ".grid.vertical { stroke: #e5e5e5; }",
        ".axis { stroke: #333; stroke-width: 1.2; }",
        ".series { stroke-width: 2.2; stroke-linecap: round; stroke-linejoin: round; }",
        ".reference { stroke: #555; stroke-width: 1.8; stroke-dasharray: 7 5; }",
        "</style>",
        "</defs>",
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="white"/>',
    ]
    svg.extend(draw_panel(spatial, series, "spatial", PANEL_X[0], "spatial"))
    svg.extend(draw_panel(temporal, series, "temporal", PANEL_X[1], "temporal"))
    svg.append("</svg>")
    output_path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def generate_figures(
    results_path,
    output_dir,
    velocity_name="velocity.svg",
    pressure_name="pressure.svg",
):
    rows = load_results(results_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_figure(
        rows,
        output_dir / velocity_name,
        [("L2_ux", "u_x"), ("L2_uy", "u_y"), ("L2_uz", "u_z")],
    )
    write_figure(
        rows,
        output_dir / pressure_name,
        [("L2_p", "pressure")],
    )
    print(f"Figures written to {output_dir}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate static convergence figures from a results CSV."
    )
    parser.add_argument(
        "results",
        nargs="?",
        type=Path,
        default=RESULTS,
        help="input CSV (default: build/convergence/results.csv)",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="figure directory (default: docs/convergence)",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="prefix added to the output filenames",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    generate_figures(
        args.results,
        args.output_dir,
        velocity_name=f"{args.prefix}velocity.svg",
        pressure_name=f"{args.prefix}pressure.svg",
    )


if __name__ == "__main__":
    main()
