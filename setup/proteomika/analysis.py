from typing import Dict, Set
import numpy as np
import pandas as pd


Comparisons = {
    "M_MUT_vs_WT": (
        "Log2FC_M_WT_M_MUT",
        "negLog10pVal_M_WT_M_MUT",
        "Samce: MUT vs WT",
    ),
    "F_MUT_vs_WT": (
        "Log2FC_F_WT_F_MUT",
        "negLog10pVal_F_WT_F_MUT",
        "Samice: MUT vs WT",
    ),
    "M_WT_vs_F_WT": (
        "Log2FC_M_WT_F_WT",
        "negLog10pVal_M_WT_F_WT",
        "Samce WT vs Samice WT",
    ),
    "M_MUT_vs_F_MUT": (
        "Log2FC_M_MUT_F_MUT",
        "negLog10pVal_M_MUT_F_MUT",
        "Samce MUT vs Samice MUT",
    ),
}


def classify_comparison(
    df: pd.DataFrame,
    comparison: str,
    fc_threshold: float = 1.2,
    p_threshold: float = 0.05,
) -> pd.DataFrame:

    if comparison not in Comparisons:
        raise ValueError(
            f"Nieznane porównanie: '{comparison}'. "
            f"Dostępne: {list(Comparisons)}"
        )

    fc_col, p_col, _ = Comparisons[comparison]
    neg_log10_threshold = -np.log10(p_threshold)

    result = pd.DataFrame({
        "Gene names":  df["Gene names"],
        "Log2FC":      df[fc_col].astype(float),
        "neg_log10_p": df[p_col].astype(float),
    }).dropna(subset=["Log2FC", "neg_log10_p"])

    sig = result["neg_log10_p"] > neg_log10_threshold
    result["status"] = np.where(
        sig & (result["Log2FC"] >  fc_threshold), "up",
        np.where(
            sig & (result["Log2FC"] < -fc_threshold), "down",
            "not_significant",
        ),
    )

    is_sig = result["status"].isin(("up", "down"))
    result["label"] = result["Gene names"].where(is_sig, "")

    return result.reset_index(drop=True)


def get_significant_sets(
    df: pd.DataFrame,
    fc_threshold: float = 1.2,
    p_threshold: float = 0.05,
) -> Dict[str, Dict[str, Set[str]]]:


    result = {}
    for comp in Comparisons:
        classified = classify_comparison(df, comp, fc_threshold, p_threshold)
        result[comp] = {
            "up":   set(classified.loc[classified["status"] == "up",   "Gene names"]),
            "down": set(classified.loc[classified["status"] == "down",  "Gene names"]),
            "all":  set(classified.loc[classified["status"] != "not_significant", "Gene names"]),
        }
    return result
