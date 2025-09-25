#!/usr/bin/env bash
set -euo pipefail

# Generate dynamic solve input files for combinations of peak_val, shear_traction, and Dc.
# - Adjust arrays below to change the parameter grid.
# - Uses the base file in this directory and writes outputs here.

BASE_FILE="dynamic_solve_alpha0d5_dsigma10_Dc0d8_db_base.i"

# Fixed residual strength used to compute dsigma (MPa offset)
RESIDUAL_SIGMA="60e6"

# Parameter grids
PEAK_VALS=("0" "0.3" "0.5")
SHEAR_TRACTIONS=("63e6" "67e6" "70e6" "73e6" "76e6")
DCS=("0.4" "0.8" "1.2")

# --- Helpers ---
script_dir() { cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd; }

format_decimal_for_name() {
  local v="$1"
  if [[ "$v" == *.* ]]; then
    echo "${v/./d}"
  else
    echo "$v"
  fi
}

to_mpa_number() {
  local s="$1"
  echo "${s%e6}"
}

main() {
  local dir; dir="$(script_dir)"
  cd "$dir"

  if [[ ! -f "$BASE_FILE" ]]; then
    echo "Base file '$BASE_FILE' not found in $dir" >&2
    exit 1
  fi

  local resid_mpa; resid_mpa="$(to_mpa_number "$RESIDUAL_SIGMA")"

  echo "Generating dynamic cases from base: $BASE_FILE"

  for peak in "${PEAK_VALS[@]}"; do
    local alpha_str; alpha_str="$(format_decimal_for_name "$peak")"
    for tau in "${SHEAR_TRACTIONS[@]}"; do
      local tau_mpa; tau_mpa="$(to_mpa_number "$tau")"
      local dsigma; dsigma=$(( tau_mpa - resid_mpa ))
      local static_mesh_rel="../static_code_files/static_solve_alpha${alpha_str}_dsigma${dsigma}_out.e"

      for dc in "${DCS[@]}"; do
        local dc_str; dc_str="$(format_decimal_for_name "$dc")"
        local out_file="dynamic_solve_alpha${alpha_str}_dsigma${dsigma}_Dc${dc_str}_db.i"

        # Perform substitutions:
        # - peak_val = <peak>
        # - Dc = <dc>
        # - mesh = '../static_code_files/static_solve_alpha{alpha}_dsigma{dsigma}_out.e'
        sed -E \
          -e "s/^([[:space:]]*peak_val[[:space:]]*=[[:space:]]*).*/\\1${peak}/" \
          -e "s/^([[:space:]]*Dc[[:space:]]*=[[:space:]]*).*/\\1${dc} #characteristic length (m)/" \
          -e "s|^([[:space:]]*mesh[[:space:]]*=[[:space:]]*).*$|\\1'${static_mesh_rel}'|" \
          "$BASE_FILE" > "$out_file"

        echo "  Wrote: $out_file (peak_val=$peak, shear_traction=$tau, Dc=$dc)"
      done
    done
  done

  echo "Done."
}

main "$@"

