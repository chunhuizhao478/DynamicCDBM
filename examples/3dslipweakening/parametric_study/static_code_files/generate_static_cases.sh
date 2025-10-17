#!/usr/bin/env bash
set -euo pipefail

# Generate static solve input files for combinations of peak_val and shear_traction.
# - Adjust the arrays below to change cases.
# - Files are written into this directory.

# Base template to copy/modify (must exist in this directory). Use the basefile
# so it is never overwritten by generation.
TEMPLATE_FILE="static_solve_alpha0_dsigma10_basefile.i"

# Residual strength (fixed for all cases)
RESIDUAL_SIGMA="60e6"

# Parameter grids (edit as needed)
# PEAK_VALS=("0" "0.3" "0.5")
PEAK_VALS=("0.15")
SHEAR_TRACTIONS=("63e6" "67e6" "70e6" "73e6" "76e6")

# --- Helpers ---
script_dir() { cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd; }

# Format alpha part for filename, e.g. 0 -> 0, 0.3 -> 0d3
format_alpha_for_name() {
  local v="$1"
  if [[ "$v" == *.* ]]; then
    echo "${v/./d}"
  else
    echo "$v"
  fi
}

# Convert strings like 70e6 -> 70 (in MPa units)
to_mpa_number() {
  local s="$1"
  echo "${s%e6}"
}

main() {
  local dir; dir="$(script_dir)"
  cd "$dir"

  if [[ ! -f "$TEMPLATE_FILE" ]]; then
    echo "Template file '$TEMPLATE_FILE' not found in $dir" >&2
    exit 1
  fi

  local resid_mpa; resid_mpa="$(to_mpa_number "$RESIDUAL_SIGMA")"

  echo "Generating cases from template: $TEMPLATE_FILE"

  for peak in "${PEAK_VALS[@]}"; do
    local alpha_str; alpha_str="$(format_alpha_for_name "$peak")"
    for tau in "${SHEAR_TRACTIONS[@]}"; do
      local tau_mpa; tau_mpa="$(to_mpa_number "$tau")"
      # dsigma in MPa: (tau - residual)/1e6; since inputs are in e6, subtract integers
      local dsigma; dsigma=$(( tau_mpa - resid_mpa ))

      local out_file="static_solve_alpha${alpha_str}_dsigma${dsigma}.i"

      # Create the file by substituting the parameters in the correct locations
      # - peak_val line
      # - shear_traction line (keep the comment format)
      # Only modify the global parameters prior to the first [Mesh] block.
      # BSD sed compatibility: apply subs on lines NOT in the range [Mesh]..EOF.
      sed -E \
        -e "/^\\[Mesh\\]/,99999! s/^([[:space:]]*peak_val[[:space:]]*=[[:space:]]*).*/\\1${peak}/" \
        -e "/^\\[Mesh\\]/,99999! s/^([[:space:]]*shear_traction[[:space:]]*=[[:space:]]*).*/\\1${tau} #Pa, shear traction/" \
        "$TEMPLATE_FILE" > "$out_file"

      echo "  Wrote: $out_file (peak_val=$peak, shear_traction=$tau)"
    done
  done

  echo "Done."
}

main "$@"
