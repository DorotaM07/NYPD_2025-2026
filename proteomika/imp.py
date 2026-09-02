import pandas as pd


Required = [
    "Gene names",
    "Log2FC_M_WT_M_MUT",
    "negLog10pVal_M_WT_M_MUT",
    "Log2FC_F_WT_F_MUT",
    "negLog10pVal_F_WT_F_MUT",
]


def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, decimal=",")
    missing = [c for c in Required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Brak wymaganych kolumn: {missing}.\n"
            f"Dostępne kolumny: {list(df.columns)}"
        )
    return df
