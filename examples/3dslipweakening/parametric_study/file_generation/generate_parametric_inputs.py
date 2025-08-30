#!/usr/bin/env python3
"""
Generate parametric MOOSE input files by varying any parameters you want.

Built-ins (backwards compatible):
- --peaks to sweep nucleation peak_shear_value (strings like 81e6)
- --dcs to sweep Slip-Weakening Dc (floats like 0.4 0.6 0.8)

General and easier way (recommended):
- --sweep 'block:key=v1,v2,v3'   to sweep a parameter in a named block
- --set   'block:key=value'      to set a constant override for all cases
- --block 'name=header|footer'   to register a custom block using header/footer regexes
    Blocks provided by default:
        slip_weakening = '##Slip weakening parameters##' ... '##-------------------------##'
        nucleation     = '#nucleation parameters' ... '##------------------------------------------------------------------##'

Filename control:
- By default, filenames append tokens for all swept keys, e.g. root_peak81e6_Dc0p8.i
- You can customize with --filename-pattern, e.g. '{root}_Dc{Dc}_peak{peak_shear_value}{ext}'

Examples:
    # Backwards compatible (peaks x dcs cartesian product)
    python generate_parametric_inputs.py \
        --base-file ../alpha/dynamic_solve_alpha0_dsigma0_elastic.i \
        --output-dir ./generated_files \
        --peaks 81e6 83e6 85e6 \
        --dcs 0.4 0.6 0.8

    # Generic: sweep and set any keys in known blocks
    python generate_parametric_inputs.py \
        --sweep 'nucleation:peak_shear_value=81e6,83e6,85e6' \
        --sweep 'slip_weakening:Dc=0.4,0.6,0.8' \
        --set   'slip_weakening:mu_s=0.65' \
        --filename-pattern '{root}_peak{peak_shear_value}_Dc{Dc}{ext}'

Notes:
- Values are treated as literals and inserted as-is (strings). For floats, compact formatting is used when necessary.
- If you need another block, register it with --block 'name=header|footer'.
"""
from __future__ import annotations

import argparse
import itertools
import os
import re
import sys
from typing import Dict, Iterable, List, Tuple


SLIP_WEAKENING_HEADER = r"##Slip weakening parameters##"
SLIP_WEAKENING_FOOTER = r"##-------------------------##"
NUCLEATION_HEADER = r"#nucleation parameters"
NUCLEATION_FOOTER = r"##------------------------------------------------------------------##"

# Default block aliases
DEFAULT_BLOCKS: Dict[str, Tuple[str, str]] = {
    "slip_weakening": (SLIP_WEAKENING_HEADER, SLIP_WEAKENING_FOOTER),
    "nucleation": (NUCLEATION_HEADER, NUCLEATION_FOOTER),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate parametric MOOSE input files.")
    parser.add_argument(
        "--base-file",
        default="../alpha/dynamic_solve_alpha0_dsigma0_elastic.i",
        help="Path to the base .i file to read.",
    )
    parser.add_argument(
        "--output-dir",
        default="./generated_files",
        help="Directory to write generated .i files.",
    )
    # Backwards compatible shortcuts
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
    # Flexible interfaces
    parser.add_argument(
        "--sweep",
        dest="sweeps",
        action="append",
        default=[],
        help=(
            "Sweep definition. Format: 'block:key=v1,v2,...'. "
            "Use block names: nucleation, slip_weakening, or define custom via --block."
        ),
    )
    parser.add_argument(
        "--set",
        dest="sets",
        action="append",
        default=[],
        help=(
            "Constant override. Format: 'block:key=value'. Applied to all generated files."
        ),
    )
    parser.add_argument(
        "--block",
        dest="custom_blocks",
        action="append",
        default=[],
        help=(
            "Register a custom block alias. Format: 'name=header|footer' (regex, anchored per line)."
        ),
    )
    parser.add_argument(
        "--filename-pattern",
        default=None,
        help=(
            "Filename pattern using {root}, {ext}, and any parameter keys as fields. "
            "Example: '{root}_Dc{Dc}_peak{peak_shear_value}{ext}'"
        ),
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

def replace_assignments_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<prefix>\s*{re.escape(key)}\s*=\s*)(?P<rhs>[^#\n]+?)(?P<suffix>\s*(#.*)?)$",
        flags=re.MULTILINE,
    )
    new_text, count = pattern.subn(rf"\g<prefix>{value}\g<suffix>", text)
    if count == 0:
        raise ValueError(f"Assignment for '{key}' not found in file (global search).")
    return new_text


def apply_param_updates(text: str, params_by_block: Dict[str, Dict[str, str]], blocks: Dict[str, Tuple[str, str]]) -> str:
    updated = text
    # First handle named blocks
    for block_name, kvs in params_by_block.items():
        if block_name == "global":
            # Defer global; do later
            continue
        if block_name not in blocks:
            raise ValueError(f"Unknown block alias '{block_name}'. Provide it via --block.")
        header, footer = blocks[block_name]
        start, end = find_block(updated, header, footer)
        block_text = updated[start:end]
        new_block = block_text
        for k, v in kvs.items():
            new_block = replace_assignment_in_block(new_block, k, v)
        updated = updated[:start] + new_block + updated[end:]
    # Then handle any global replacements
    global_kvs = params_by_block.get("global", {})
    for k, v in global_kvs.items():
        updated = replace_assignments_global(updated, k, v)
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


def render_filename(base_path: str, params: Dict[str, str], pattern: str | None) -> str:
    base_dir, base_name = os.path.dirname(base_path), os.path.basename(base_path)
    root, ext = os.path.splitext(base_name)
    if pattern:
        # Provide sanitized tokens for all keys
        fmt_params = {k: sanitize_for_filename(v) for k, v in params.items()}
        fmt_params.update({"root": root, "ext": ext})
        new_name = pattern.format(**fmt_params)
    else:
        # Default pattern: append tokens for sorted keys
        parts = [root]
        for k in sorted(params.keys()):
            parts.append(f"{k}{sanitize_for_filename(params[k])}")
        new_name = "_".join(parts) + ext
    return os.path.join(base_dir, new_name)


def parse_keyval_list(arg: str) -> Tuple[str, str, List[str]]:
    """Parse 'block:key=v1,v2,...' into (block, key, [values]).
    If 'block' is omitted (':key=...'), treat as global.
    """
    if ":" in arg.split("=", 1)[0]:
        left, values_str = arg.split("=", 1)
        block, key = left.split(":", 1)
        block = block.strip() or "global"
        key = key.strip()
    else:
        # Allow 'key=v1,v2' shorthand -> global block
        key, values_str = arg.split("=", 1)
        block = "global"
        key = key.strip()
    values = [v.strip() for v in values_str.split(",") if v.strip() != ""]
    if not key or not values:
        raise ValueError(f"Invalid parameter specification: {arg}")
    return block, key, values


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

    # Build blocks mapping (defaults + custom)
    blocks = dict(DEFAULT_BLOCKS)
    for b in args.custom_blocks:
        try:
            name, rest = b.split("=", 1)
            header, footer = rest.split("|", 1)
        except ValueError:
            print(f"ERROR: Invalid --block format: {b}. Expected 'name=header|footer'", file=sys.stderr)
            return 2
        blocks[name.strip()] = (header.strip(), footer.strip())

    # Collect sweeps and constant sets
    sweeps: Dict[str, Dict[str, List[str]]] = {}
    sets: Dict[str, Dict[str, str]] = {}

    # Backwards compatibility: map --peaks and --dcs into sweeps if provided
    if peaks:
        sweeps.setdefault("nucleation", {})["peak_shear_value"] = list(peaks)
    if dcs:
        # Render as compact strings
        sweeps.setdefault("slip_weakening", {})["Dc"] = [format_number_literal(x) for x in dcs]

    # Parse generic --sweep flags
    try:
        for s in args.sweeps:
            block, key, values = parse_keyval_list(s)
            sweeps.setdefault(block, {})[key] = values
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    # Parse generic --set flags
    try:
        for s in args.sets:
            block, key, values = parse_keyval_list(s)
            if len(values) != 1:
                print(f"ERROR: --set requires a single value, got: {s}", file=sys.stderr)
                return 2
            sets.setdefault(block, {})[key] = values[0]
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    # Build the cartesian product across all swept keys
    # Flatten sweeps into a list of (block, key, [values])
    sweep_items: List[Tuple[str, str, List[str]]] = []
    for block, kvs in sweeps.items():
        for key, values in kvs.items():
            sweep_items.append((block, key, values))

    if not sweep_items:
        print("ERROR: No sweep parameters specified. Use --sweep/--peaks/--dcs.", file=sys.stderr)
        return 2

    # Prepare product
    value_lists = [vals for (_, _, vals) in sweep_items]
    products = list(itertools.product(*value_lists))

    made = []
    for combo in products:
        # Construct params_by_block for this combination, starting with constants
        params_by_block: Dict[str, Dict[str, str]] = {blk: kv.copy() for blk, kv in sets.items()}
        flat_params: Dict[str, str] = {k: v for blk, kv in sets.items() for k, v in kv.items()}

        for (block, key, _), value in zip(sweep_items, combo):
            params_by_block.setdefault(block, {})[key] = value
            flat_params[key] = value

        try:
            new_text = apply_param_updates(base_text, params_by_block, blocks)
        except ValueError as e:
            print(f"ERROR updating parameters for {flat_params}: {e}", file=sys.stderr)
            return 3

        # Name and path
        out_name = os.path.basename(render_filename(base_file, flat_params, args.filename_pattern))
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
