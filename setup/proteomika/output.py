import cProfile
import io
import os
import pstats
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go


def save_volcano(fig: go.Figure, output_dir: str, name: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"volcano_{name}.html")
    fig.write_html(path, include_plotlyjs="cdn")
    return path


def save_venn(fig: plt.Figure, output_dir: str, name: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"venn_{name}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def save_significant_csv(
    classified: pd.DataFrame,
    output_dir: str,
    name: str,
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    sig = classified[classified["status"].isin(["up", "down"])].copy()
    path = os.path.join(output_dir, f"significant_{name}.csv")
    sig.to_csv(path, index=False)
    return path


def save_report(lines: list, output_dir: str, filename: str = "report.txt") -> str:
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


def save_profile_report(profiler: cProfile.Profile, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    buf = io.StringIO()
    ps = pstats.Stats(profiler, stream=buf).sort_stats("cumulative")
    ps.print_stats(30)

    path = os.path.join(output_dir, "profile_report.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Raport profilera\n")
        fh.write("Sortowanie: cumulative\n\n")
        fh.write(buf.getvalue())
    return path


def print_summary(lines: list) -> None:
    for line in lines:
        print(line)
