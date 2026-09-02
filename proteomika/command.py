import argparse
import cProfile
import os
import sys

from proteomika.imp import load_data
from proteomika.analysis import Comparisons, classify_comparison, get_significant_sets
from proteomika.volcano_plot import make_volcano
from proteomika.Venn_plot import make_venn, venn_summary_text
from proteomika.output import (
    save_volcano,
    save_venn,
    save_significant_csv,
    save_report,
    save_profile_report,
    print_summary,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="proteomika",
        description="Analiza proteomiki: volcano ploty i diagramy Venna.",
    )
    parser.add_argument(
        "--input",
        required=True,
        metavar="Plik.csv",
        help="Plik CSV z danymi proteomicznymi.",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="Katalog",
        help="Katalog wyjsciowy na pliki HTML, PNG i CSV.",
    )
    parser.add_argument(
        "--fc",
        type=float,
        default=1.2,
        metavar="N",
        help="Prog log2FC (domyslnie: 1.2).",
    )
    parser.add_argument(
        "--pvalue",
        type=float,
        default=0.05,
        metavar="N",
        help="Prog p-value (domyslnie: 0.05).",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Uruchamianie z profilerem.",
    )
    return parser


def run(args: argparse.Namespace) -> None:
    report_lines = [
        "Raport analizy proteomicznej",
        f"Plik wejsciowy: {args.input}",
        f"Prog log2FC:    {args.fc}",
        f"Prog p-value:   {args.pvalue}",
        "",
    ]

    print(f"Wczytuje dane: {args.input}")
    try:
        df = load_data(args.input)
    except (ValueError, FileNotFoundError) as exc:
        print(f"błąd: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"  Wczytano {len(df)} bialek.")

    print("\nTworzenie volcano plotow...")
    for comp_key, (_, _, label) in Comparisons.items():
        classified = classify_comparison(df, comp_key, args.fc, args.pvalue)

        up_n   = (classified["status"] == "up").sum()
        down_n = (classified["status"] == "down").sum()
        report_lines.append(f"{label}:")
        report_lines.append(f"  Up-regulated:   {up_n}")
        report_lines.append(f"  Down-regulated: {down_n}")
        report_lines.append(f"  Istotne: {up_n + down_n}")
        report_lines.append("")

        fig = make_volcano(classified, args.fc, args.pvalue, title=label)
        vpath = save_volcano(fig, args.output, comp_key)
        print(f"  [{comp_key}] -> {vpath}")

        cpath = save_significant_csv(classified, args.output, comp_key)
        print(f"  [{comp_key}] -> {cpath}")

    print("\nTworzenie wykresów Venna...")
    sig_sets = get_significant_sets(df, args.fc, args.pvalue)

    for direction in ("up", "down"):
        venn_fig, info = make_venn(sig_sets, direction)
        venn_path = save_venn(venn_fig, args.output, f"M_vs_F_{direction}")
        print(f"  [Venn {direction}] -> {venn_path}")

        summary = venn_summary_text(info, direction)
        report_lines.append(summary)

    rpath = save_report(report_lines, args.output)
    print(f"\nRaport: {rpath}")
    print_summary(report_lines)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not args.profile:
        run(args)
        return

    print("Profilowanie...")
    profiler = cProfile.Profile()
    profiler.enable()
    run(args)
    profiler.disable()

    ppath = save_profile_report(profiler, args.output)
    print(f"Raport profilera: {ppath}")


if __name__ == "__main__":
    main()