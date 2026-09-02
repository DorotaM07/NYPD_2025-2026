import io
import os
import numpy as np
import pandas as pd
import pytest

from proteomika.imp import load_data, Required
from proteomika.analysis import classify_comparison, get_significant_sets, Comparisons
from proteomika.volcano_plot import make_volcano
from proteomika.Venn_plot import make_venn, venn_summary_text
from proteomika.output import save_volcano, save_venn, save_significant_csv, save_report

@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "Gene names":              ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"],
        "Log2FC_M_WT_M_MUT":      [ 2.0,    -2.0,     0.1,    -0.1,     1.5],
        "negLog10pVal_M_WT_M_MUT":[ 2.0,     2.0,     0.5,     0.5,     2.0],
        "Log2FC_F_WT_F_MUT":      [ 2.0,     0.1,    -2.0,    -0.1,     0.1],
        "negLog10pVal_F_WT_F_MUT":[ 2.0,     0.5,     2.0,     0.5,     0.5],
        "Log2FC_M_WT_F_WT":       [ 0.1,     0.1,     0.1,     0.1,     0.1],
        "negLog10pVal_M_WT_F_WT": [ 0.5,     0.5,     0.5,     0.5,     0.5],
        "Log2FC_M_MUT_F_MUT":     [ 0.1,     0.1,     0.1,     0.1,     0.1],
        "negLog10pVal_M_MUT_F_MUT":[ 0.5,    0.5,     0.5,     0.5,     0.5],
    })


@pytest.fixture
def sample_csv(tmp_path, sample_df):
    path = tmp_path / "test.csv"
    sample_df.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def classified_m(sample_df):
    return classify_comparison(sample_df, "M_MUT_vs_WT", fc_threshold=1.0, p_threshold=0.05)

class TestLoadData:
    def test_loads_csv(self, sample_csv):
        df = load_data(sample_csv)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 5

    def test_missing_column_raises(self, tmp_path):
        bad = pd.DataFrame({"Gene names": ["A"], "Log2FC_M_WT_M_MUT": [1.0]})
        path = tmp_path / "bad.csv"
        bad.to_csv(path, index=False)
        with pytest.raises(ValueError, match="Brak wymaganych kolumn"):
            load_data(str(path))

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            load_data("nieistniejacy_plik.csv")

    def test_required_columns_present(self, sample_csv):
        df = load_data(sample_csv)
        for col in Required:
            assert col in df.columns


class TestClassifyComparison:
    def test_returns_dataframe(self, sample_df):
        result = classify_comparison(sample_df, "M_MUT_vs_WT")
        assert isinstance(result, pd.DataFrame)

    def test_required_output_columns(self, classified_m):
        for col in ("Gene names", "Log2FC", "neg_log10_p", "status", "label"):
            assert col in classified_m.columns

    def test_up_count(self, classified_m):
        assert (classified_m["status"] == "up").sum() == 2

    def test_down_count(self, classified_m):
        assert (classified_m["status"] == "down").sum() == 1

    def test_label_empty_for_not_significant(self, classified_m):
        not_sig = classified_m[classified_m["status"] == "not_significant"]
        assert (not_sig["label"] == "").all()

    def test_label_filled_for_significant(self, classified_m):
        sig = classified_m[classified_m["status"].isin(["up", "down"])]
        assert (sig["label"] != "").all()

    def test_unknown_comparison_raises(self, sample_df):
        with pytest.raises(ValueError, match="Błąd- nieznane porównanie"):
            classify_comparison(sample_df, "Nieznane")

    def test_fc_threshold_respected(self, sample_df):
        result = classify_comparison(sample_df, "M_MUT_vs_WT", fc_threshold=3.0)
        assert (result["status"] == "not_significant").all()

    def test_does_not_modify_input(self, sample_df):
        cols_before = set(sample_df.columns)
        classify_comparison(sample_df, "M_MUT_vs_WT")
        assert set(sample_df.columns) == cols_before

class TestGetSignificantSets:
    def test_returns_all_comparisons(self, sample_df):
        sets = get_significant_sets(sample_df)
        assert set(sets.keys()) == set(Comparisons.keys())

    def test_each_comparison_has_up_down_all(self, sample_df):
        sets = get_significant_sets(sample_df)
        for comp in sets.values():
            assert "up" in comp and "down" in comp and "all" in comp

    def test_up_set_correct_M(self, sample_df):
        sets = get_significant_sets(sample_df)
        assert "GeneA" in sets["M_MUT_vs_WT"]["up"]
        assert "GeneE" in sets["M_MUT_vs_WT"]["up"]

    def test_down_set_correct_M(self, sample_df):
        sets = get_significant_sets(sample_df)
        assert "GeneB" in sets["M_MUT_vs_WT"]["down"]

    def test_all_is_union_of_up_and_down(self, sample_df):
        sets = get_significant_sets(sample_df)
        for comp in sets.values():
            assert comp["all"] == comp["up"] | comp["down"]

class TestMakeVolcano:
    def test_returns_figure(self, classified_m):
        import plotly.graph_objects as go
        fig = make_volcano(classified_m, 1.0, 0.05)
        assert isinstance(fig, go.Figure)

    def test_has_traces(self, classified_m):
        fig = make_volcano(classified_m, 1.0, 0.05)
        assert len(fig.data) > 0

    def test_custom_title(self, classified_m):
        fig = make_volcano(classified_m, 1.0, 0.05, title="Test")
        assert fig.layout.title.text == "Test"

class TestMakeVenn:
    def test_returns_figure_and_info(self, sample_df):
        sets = get_significant_sets(sample_df)
        fig, info = make_venn(sets, "up")
        import matplotlib.pyplot as plt
        assert isinstance(fig, plt.Figure)
        assert "only_M" in info and "only_F" in info and "shared" in info

    def test_up_sets_correct(self, sample_df):
        sets = get_significant_sets(sample_df)
        _, info = make_venn(sets, "up")
        assert "GeneA" in info["shared"]
        assert "GeneE" in info["only_M"]

    def test_down_sets_correct(self, sample_df):
        sets = get_significant_sets(sample_df)
        _, info = make_venn(sets, "down")
        assert "GeneB" in info["only_M"]
        assert "GeneC" in info["only_F"]

    def test_invalid_direction_raises(self, sample_df):
        sets = get_significant_sets(sample_df)
        with pytest.raises(ValueError, match="direction"):
            make_venn(sets, "invalid")

    def test_venn_summary_text_contains_direction(self, sample_df):
        sets = get_significant_sets(sample_df)
        _, info = make_venn(sets, "up")
        text = venn_summary_text(info, "up")
        assert "Up-regulated" in text

class TestOutput:
    def test_save_volcano_creates_html(self, classified_m, tmp_path):
        fig = make_volcano(classified_m, 1.0, 0.05)
        path = save_volcano(fig, str(tmp_path), "test")
        assert os.path.exists(path)
        assert path.endswith(".html")

    def test_save_venn_creates_png(self, sample_df, tmp_path):
        sets = get_significant_sets(sample_df)
        fig, _ = make_venn(sets, "up")
        path = save_venn(fig, str(tmp_path), "test_up")
        assert os.path.exists(path)
        assert path.endswith(".png")

    def test_save_significant_csv_correct_count(self, classified_m, tmp_path):
        path = save_significant_csv(classified_m, str(tmp_path), "test")
        result = pd.read_csv(path)
        assert len(result) == 3

    def test_save_report_creates_file(self, tmp_path):
        path = save_report(["linia 1", "linia 2"], str(tmp_path))
        assert os.path.exists(path)
        content = open(path).read()
        assert "linia 1" in content
