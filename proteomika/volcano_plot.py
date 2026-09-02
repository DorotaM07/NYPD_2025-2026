import numpy as np
import pandas as pd
import plotly.graph_objects as go


Colors = {
    "not_significant": "#aaaaaa",
    "up":              "#d62728",
    "down":            "#1f77b4",
}

Groups = [
    ("Not significant", "not_significant", "circle"),
    ("Up-regulated",    "up",              "circle"),
    ("Down-regulated",  "down",            "circle"),
]


def make_volcano(
    classified: pd.DataFrame,
    fc_threshold: float,
    p_threshold: float,
    title: str = "Volcano Plot",
) -> go.Figure:
    
    neg_log10_thr = -np.log10(p_threshold)
    fig = go.Figure()

    for group_label, status, symbol in Groups:
        subset = classified[classified["status"] == status]
        if subset.empty:
            continue

        fig.add_trace(go.Scatter(
            x=subset["Log2FC"],
            y=subset["neg_log10_p"],
            mode="markers+text",
            name=group_label,
            text=subset["label"].fillna(""),
            textposition="top center",
            textfont=dict(size=9),
            marker=dict(
                color=Colors[status],
                size=7,
                symbol=symbol,
                opacity=0.8,
            ),
            customdata=subset[["Gene names", "Log2FC", "neg_log10_p"]].values,
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "Log2FC: %{customdata[1]:.3f}<br>"
                "-log10(p): %{customdata[2]:.3f}<br>"
                "<extra></extra>"
            ),
        ))

    fig.add_hline(
        y=neg_log10_thr,
        line_dash="dash",
        line_color="black",
        opacity=0.5,
        annotation_text=f"p={p_threshold}",
        annotation_position="right",
    )
    fig.add_vline(x= fc_threshold, line_dash="dash", line_color="black", opacity=0.5)
    fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="black", opacity=0.5)

    fig.update_layout(
        title=title,
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10(p-value)",
        legend_title="Status",
        template="plotly_white",
        height=650,
    )
    return fig
