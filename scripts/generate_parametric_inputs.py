#!/usr/bin/env python3
"""
Generate parametric MOOSE input files by varying peak_shear_value and Dc.

- Reads a base input file (MOOSE .i) that contains the following sections:
  ##Slip weakening parameters## ... Dc = <value> ... ##-------------------------##
  #nucleation parameters ... peak_shear_value = <value> ... ##------------------------------------------------------------------##

- Writes one derived input file for each combination of provided peaks and Dcs.
- New filenames embed the values, e.g.,
  dynamic_solve_alpha0_dsigma0_elastic_peak81e6_Dc0p8.i

Usage examples:
  python scripts/generate_parametric_inputs.py \
    --base-file examples/3dslipweakening/parametric_study/alpha/dynamic_solve_alpha0_dsigma0_elastic.i \
    --output-dir examples/3dslipweakening/parametric_study/alpha/generated \
    --peaks 81e6 83e6 85e6 \
    --dcs 0.4 0.6 0.8

Notes:
- Peaks are treated as strings to preserve scientific notation like 81e6.
- Dc values are floats but rendered compactly; decimals are "sanitized" to filenames by replacing '.' with 'p'.
"""
from __future__ import annotations

import argparse
import itertools
import os
import re
import sys
from typing import Iterable, Tuple


SLIP_WEAKENING_HEADER = r"##Slip weakening parameters##"
SLIP_WEAKENING_FOOTER = r"##-------------------------##"
NUCLEATION_HEADER = r"#nucleation parameters"
NUCLEATION_FOOTER = r"##------------------------------------------------------------------##"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate parametric MOOSE input files.")
    parser.add_argument(
        "--base-file",
        default="examples/3dslipweakening/parametric_study/alpha/dynamic_solve_alpha0_dsigma0_elastic.i",
        help="Path to the base .i file to read.",
    )
    parser.add_argument(
        "--output-dir",
        default="examples/3dslipweakening/parametric_study/alpha/generated",
        help="Directory to write generated .i files.",
    )
    parser.add_argument(
        "--peaks",
        nargs="+",
        type=str,
        required=False,
        default=["81e6", "83e6", "85e6"],
        help="List of peak_shear_value values (strings, e.g., 81e6 83e6 85e6).",
    )
    parser.add_argument(
        "--dcs",
        nargs="+",
        type=float,
        required=False,
        default=[0.4, 0.6, 0.8],
        help="List of Dc values (floats, e.g., 0.4 0.6 0.8).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse and show what would be generated without writing files.",
    )
    return parser.parse_args()


def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def write_text(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def find_block(text: str, header_pat: str, footer_pat: str) -> Tuple[int, int]:
    """Return (start_index, end_index) spanning the block including header line up to just before footer line.
    Raises ValueError if not found.
    """
    header_match = re.search(rf"^\s*{header_pat}\s*$", text, flags=re.MULTILINE)
    if not header_match:
        raise ValueError(f"Header not found: {header_pat}")
    footer_match = re.search(rf"^\s*{footer_pat}\s*$", text[header_match.end():], flags=re.MULTILINE)
    if not footer_match:
        raise ValueError(f"Footer not found after header {header_pat}: {footer_pat}")
    start = header_match.start()
    end = header_match.end() + footer_match.start()
    return start, end


def replace_assignment_in_block(block_text: str, key: str, value: str) -> str:
    """Replace a single-line assignment of the form 'key = <anything up to # or EOL>' inside block_text.
    Preserves any trailing inline comment.
    Raises ValueError if the key assignment isn't found exactly once.
    """
    pattern = re.compile(
        rf"^(?P<prefix>\s*{re.escape(key)}\s*=\s*)(?P<rhs>[^#\n]+?)(?P<suffix>\s*(#.*)?)$",
        flags=re.MULTILINE,
    )
    new_block, count = pattern.subn(rf"\g<prefix>{value}\g<suffix>", block_text)
    if count == 0:
        raise ValueError(f"Assignment for '{key}' not found in its block.")
    if count > 1:
        raise ValueError(f"Multiple assignments for '{key}' found in its block ({count}). Ambiguous.")
    return new_block


def set_parameters(text: str, peak_value: str, dc_value: float) -> str:
    # Replace Dc in the Slip weakening parameters block only
    sw_start, sw_end = find_block(text, SLIP_WEAKENING_HEADER, SLIP_WEAKENING_FOOTER)
    sw_block = text[sw_start:sw_end]
    sw_block_new = replace_assignment_in_block(sw_block, "Dc", format_number_literal(dc_value))

    # Replace peak_shear_value in the nucleation parameters block only
    nu_start, nu_end = find_block(text, NUCLEATION_HEADER, NUCLEATION_FOOTER)
    nu_block = text[nu_start:nu_end]
    nu_block_new = replace_assignment_in_block(nu_block, "peak_shear_value", peak_value)

    # Reconstruct the full text
    updated = text[:sw_start] + sw_block_new + text[sw_end:nu_start] + nu_block_new + text[nu_end:]
    return updated


def sanitize_for_filename(val: str) -> str:
    s = str(val)
    # Common sanitization: replace dot, plus, minus, spaces
    s = s.replace(".", "p").replace("+", "").replace("-", "m").replace(" ", "")
    return s


def format_number_literal(val: float) -> str:
    """Format Dc as a compact literal matching typical .i style (e.g., '0.4' not '0.400000')."""
    # Use rstrip to remove trailing zeros and dot
    s = f"{val:.10g}"
    # Ensure at least one decimal for values between 0 and 1 when appropriate (e.g., 0.4)
    # .10g already yields compact form.
    return s


def generate_names(base_path: str, peak: str, dc: float) -> str:
    base_dir, base_name = os.path.dirname(base_path), os.path.basename(base_path)
    root, ext = os.path.splitext(base_name)
    peak_token = sanitize_for_filename(peak)
    dc_token = sanitize_for_filename(format_number_literal(dc))
    new_name = f"{root}_peak{peak_token}_Dc{dc_token}{ext}"
    return os.path.join(base_dir, new_name)


def main() -> int:
    args = parse_args()

    base_file = args.base_file
    out_dir = args.output_dir
    peaks: Iterable[str] = args.peaks
    dcs: Iterable[float] = args.dcs

    if not os.path.isfile(base_file):
        print(f"ERROR: Base file not found: {base_file}", file=sys.stderr)
        return 2

    base_text = read_text(base_file)

    # We'll write outputs alongside out_dir with same filename root, but we keep name generation consistent.
    os.makedirs(out_dir, exist_ok=True)

    made = []
    for peak, dc in itertools.product(peaks, dcs):
        try:
            new_text = set_parameters(base_text, peak_value=peak, dc_value=dc)
        except ValueError as e:
            print(f"ERROR updating parameters for peak={peak}, Dc={dc}: {e}", file=sys.stderr)
            return 3

        # Name and path
        out_name = os.path.basename(generate_names(base_file, peak, dc))
        out_path = os.path.join(out_dir, out_name)

        if args.dry_run:
            print(f"Would write: {out_path}")
        else:
            write_text(out_path, new_text)
            made.append(out_path)

    if args.dry_run:
        print("Dry run complete.")
    else:
        print(f"Wrote {len(made)} files:")
        for p in made:
            print(f"  {p}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
