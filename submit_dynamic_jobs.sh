#!/usr/bin/env bash
set -euo pipefail

# Submit selected dynamic solve sbatch jobs in bulk.
# Filters by alpha, dsigma, and Dc identifiers derived from filenames
# like job_frontera_dynamic_solve_alpha0d5_dsigma10_Dc0d8_db.sbatch.

DEFAULT_JOB_SUBDIR="jobs/dynamic_solve"

usage() {
  cat <<USAGE
Usage: $(basename "$0") -a value [-a value ...] -d value [-d value ...] -c value [-c value ...] [options]

Options:
  -a value     Alpha values to include (e.g. 0, 0.3, 0d5). Repeat for multiple.
  -d value     Dsigma values to include (e.g. 3, 10, 10.0). Repeat for multiple.
  -c value     Dc values to include (e.g. 0.4, 0d8, 1.2). Repeat for multiple.
  -j dir       Directory containing generated sbatch files (default: $DEFAULT_JOB_SUBDIR relative to repo root).
  -n           Dry run; print matching sbatch commands without submitting.
  -h           Show this help message.

All three selectors (-a, -d, -c) are required. Values may be provided in either
floating form (0.5) or filename form (0d5). The script must be run from within
this repository (or with -j pointing to the correct jobs directory).
USAGE
}

script_dir() {
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd
}

# Accept values already in filename form (with 'd') or decimal form.
format_component() {
  local value="$1"
  if [[ -z "$value" ]]; then
    echo ""; return
  fi
  if [[ "$value" =~ ^[0-9]+d[0-9]+$ || "$value" =~ ^[0-9]+$ ]]; then
    echo "$value"
    return
  fi
  # Normalize leading decimal point (e.g. .5 -> 0.5)
  if [[ "$value" =~ ^\.[0-9]+$ ]]; then
    value="0${value}"
  fi
  # Trim trailing .0 sequences for integers.
  if [[ "$value" =~ ^([0-9]+)\.0+$ ]]; then
    echo "${BASH_REMATCH[1]}"
    return
  fi
  # General decimal replacement: convert '.' to 'd'.
  echo "${value//./d}"
}

format_dsigma() {
  local value="$1"
  if [[ -z "$value" ]]; then
    echo ""; return
  fi
  if [[ "$value" =~ ^[0-9]+$ ]]; then
    echo "$value"
    return
  fi
  if [[ "$value" =~ ^[0-9]+d[0-9]+$ ]]; then
    echo "$value"
    return
  fi
  if [[ "$value" =~ ^([0-9]+)\.0+$ ]]; then
    echo "${BASH_REMATCH[1]}"
    return
  fi
  if [[ "$value" =~ ^([0-9]+)\.([0-9]+)$ ]]; then
    echo "${BASH_REMATCH[1]}d${BASH_REMATCH[2]}"
    return
  fi
  echo "$value"
}

main() {
  local repo_root job_dir dry_run=false
  repo_root="$(script_dir)"
  job_dir="$repo_root/$DEFAULT_JOB_SUBDIR"

  declare -a ALPHAS=()
  declare -a DSIGMAS=()
  declare -a DCS=()

  while getopts ":a:d:c:j:hn" opt; do
    case "$opt" in
      a) ALPHAS+=("$OPTARG") ;;
      d) DSIGMAS+=("$OPTARG") ;;
      c) DCS+=("$OPTARG") ;;
      j) job_dir="$OPTARG" ;;
      n) dry_run=true ;;
      h)
        usage
        exit 0
        ;;
      :) echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
      \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    esac
  done
  shift $((OPTIND - 1))

  if [[ ${#ALPHAS[@]} -eq 0 || ${#DSIGMAS[@]} -eq 0 || ${#DCS[@]} -eq 0 ]]; then
    echo "Error: must supply at least one value for each of -a, -d, and -c." >&2
    usage
    exit 1
  fi

  # Resolve job directory (allow relative paths)
  if [[ ! "$job_dir" =~ ^/ ]]; then
    job_dir="$(cd "$job_dir" 2>/dev/null && pwd)"
  fi
  if [[ -z "$job_dir" || ! -d "$job_dir" ]]; then
    echo "Job directory not found." >&2
    exit 1
  fi

  # Normalize selectors into filename components.
  local formatted_alphas=()
  local formatted_dsigmas=()
  local formatted_dcs=()

  for alpha in "${ALPHAS[@]}"; do
    formatted_alphas+=("$(format_component "$alpha")")
  done
  for dsigma in "${DSIGMAS[@]}"; do
    formatted_dsigmas+=("$(format_dsigma "$dsigma")")
  done
  for dc in "${DCS[@]}"; do
    formatted_dcs+=("$(format_component "$dc")")
  done

  local submitted=0 missing=0 failed=0
  for alpha in "${formatted_alphas[@]}"; do
    for dsigma in "${formatted_dsigmas[@]}"; do
      for dc in "${formatted_dcs[@]}"; do
        local job_file="${job_dir}/job_frontera_dynamic_solve_alpha${alpha}_dsigma${dsigma}_Dc${dc}_db.sbatch"
        if [[ ! -f "$job_file" ]]; then
          echo "Missing sbatch file: $job_file" >&2
          ((missing++))
          continue
        fi
        if [[ "$dry_run" == true ]]; then
          echo "sbatch \"$job_file\""
          ((submitted++))
        else
          echo "Submitting $job_file"
          if sbatch "$job_file"; then
            ((submitted++))
          else
            echo "sbatch failed for $job_file" >&2
            ((failed++))
          fi
        fi
      done
    done
  done

  if [[ "$dry_run" == true ]]; then
    echo "Dry run complete. Commands listed: $submitted. Missing files: $missing." >&2
  else
    echo "Done. Submitted $submitted jobs. Missing files: $missing. sbatch failures: $failed." >&2
  fi
}

main "$@"
