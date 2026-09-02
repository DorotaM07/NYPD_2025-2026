from typing import Dict, Set, Tuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def make_venn(
    sets: Dict[str, Dict[str, Set[str]]],
    direction: str,
    title: str = "",
) -> Tuple[plt.Figure, Dict]:
    if direction not in ("up", "down"):
        raise ValueError(f"direction musi być 'up' lub 'down', otrzymano: '{direction}'")

    set_m = sets["M_MUT_vs_WT"][direction]
    set_f = sets["F_MUT_vs_WT"][direction]

    only_m = set_m - set_f
    only_f = set_f - set_m
    shared = set_m & set_f

    fig, ax = plt.subplots(figsize=(7, 5))

    venn2(
        subsets=(len(only_m), len(only_f), len(shared)),
        set_labels=("Samce MUT vs WT", "Samice MUT vs WT"),
        ax=ax,
    )

    direction_label = "up-regulated" if direction == "up" else "down-regulated"
    ax.set_title(title or f"Białka {direction_label}: Samce vs Samice MUT", fontsize=13)

    plt.tight_layout()

    info = {
        "only_M":   sorted(only_m),
        "only_F":   sorted(only_f),
        "shared":   sorted(shared),
    }
    return fig, info


def venn_summary_text(info: Dict, direction: str) -> str:
    direction_label = "Up-regulated" if direction == "up" else "Down-regulated"
    lines = [
        f"\n Venn: {direction_label} białka (M MUT vs F MUT)",
        f"Unikalne dla samców:  {len(info['only_M'])} : {info['only_M']}",
        f"Unikalne dla samic:   {len(info['only_F'])} : {info['only_F']}",
        f"Wspólne:              {len(info['shared'])} : {info['shared']}",
    ]
    return "\n".join(lines)
