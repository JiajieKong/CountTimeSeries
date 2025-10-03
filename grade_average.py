"""
A small CLI program that reads a CSV file of student grades,
calculates each student's average (treating missing values as zeros),
and writes a new CSV with columns: Name, Average.

Usage:
    python grade_average.py input.csv output.csv

Constraints:
    - Uses only the Python standard library.
    - Missing or non-numeric values are treated as 0.
"""

from __future__ import annotations

import argparse
import csv
from typing import Iterable, List, Any


def _to_float_or_zero(value: Any) -> float:
    """Convert a value to float, treating blanks or invalid as 0.0."""
    if value is None:
        return 0.0
    text = str(value).strip()
    if text == "":
        return 0.0
    try:
        return float(text)
    except (ValueError, TypeError):
        return 0.0


def calculate_average(scores: Iterable[Any]) -> float:
    """Calculate the average of scores treating missing/invalid entries as 0.0.

    The denominator is the number of provided score fields (i.e., len(scores)).
    If no scores are provided, returns 0.0.
    """
    items: List[Any] = list(scores)
    if not items:
        return 0.0
    total = sum(_to_float_or_zero(s) for s in items)
    return total / float(len(items))


def process_grades(input_csv_path: str, output_csv_path: str) -> None:
    """Read input CSV, compute per-student average, and write output CSV.

    Assumes the first column is the student's name and all remaining columns are scores.
    Missing or non-numeric values are treated as 0 and included in the denominator.
    """
    with open(input_csv_path, newline="", encoding="utf-8") as fin:
        reader = csv.reader(fin)
        try:
            header = next(reader)
        except StopIteration:
            # Empty input: write only header to output
            with open(output_csv_path, "w", newline="", encoding="utf-8") as fout:
                writer = csv.writer(fout)
                writer.writerow(["Name", "Average"])
            return

        num_score_cols = max(len(header) - 1, 0)

        with open(output_csv_path, "w", newline="", encoding="utf-8") as fout:
            writer = csv.writer(fout)
            writer.writerow(["Name", "Average"])

            for row in reader:
                # Normalize row length to at least name + num_score_cols
                if len(row) < 1 + num_score_cols:
                    row = row + [""] * (1 + num_score_cols - len(row))

                name = row[0] if row else ""
                scores = row[1:1 + num_score_cols]
                avg = calculate_average(scores)
                writer.writerow([name, avg])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Read a CSV of student grades and write Name,Average with "
            "missing values treated as zeros."
        )
    )
    parser.add_argument("input_csv", help="Path to input CSV file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    args = parser.parse_args()

    process_grades(args.input_csv, args.output_csv)


if __name__ == "__main__":
    main()

