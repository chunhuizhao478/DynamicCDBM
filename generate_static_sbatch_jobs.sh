#!/usr/bin/env bash
set -euo pipefail

# Generate Frontier sbatch job files for each static solve input deck.
# Usage: ./generate_static_sbatch_jobs.sh [-s scratch_base] [-o output_dir] [-e executable]
#  -s scratch_base : Remote scratch directory that mirrors this repo (default: /scratch1/10024/zhaochun/projects/DynamicCDBM_09202025)
#  -o output_dir   : Directory to place generated .sbatch files (default: jobs/static_solve)
#  -e executable   : Path to dynamic_cdbm executable used in the job script (default: ./dynamic_cdbm-opt)

DEFAULT_SCRATCH_BASE="/scratch1/10024/zhaochun/projects/DynamicCDBM_10022025"
DEFAULT_EXECUTABLE="./dynamic_cdbm-opt"
DEFAULT_OUTPUT_SUBDIR="jobs/static_solve"

script_dir() {
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd
}

show_help() {
  cat <<USAGE
Usage: $(basename "$0") [-s scratch_base] [-o output_dir] [-e executable]

Creates Frontier sbatch job files that mirror the structure of the template
job (job_frontera_static_solve_alpha0d5_dsigma10.sbatch) for every static
solve input deck located under:
  examples/3dslipweakening/parametric_study/static_code_files

Options:
  -s scratch_base   Remote root directory that mirrors the repository (default: $DEFAULT_SCRATCH_BASE)
  -o output_dir     Directory to receive generated sbatch files (default: jobs/static_solve)
  -e executable     Executable invoked by ibrun (default: $DEFAULT_EXECUTABLE)
  -h                Show this help and exit

Environment overrides:
  SCRATCH_BASE, OUTPUT_DIR, MOOSE_EXEC
USAGE
}

main() {
  local repo_root static_dir scratch_base output_dir executable default_output_dir
  repo_root="$(script_dir)"
  static_dir="$repo_root/examples/3dslipweakening/parametric_study/static_code_files"
  default_output_dir="$repo_root/$DEFAULT_OUTPUT_SUBDIR"

  scratch_base="${SCRATCH_BASE:-$DEFAULT_SCRATCH_BASE}"
  output_dir="${OUTPUT_DIR:-$default_output_dir}"
  executable="${MOOSE_EXEC:-$DEFAULT_EXECUTABLE}"

  local opt
  while getopts ":s:o:e:h" opt; do
    case "$opt" in
      s) scratch_base="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      e) executable="$OPTARG" ;;
      h)
        show_help
        exit 0
        ;;
      :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
      \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
  done

  if [[ ! -d "$static_dir" ]]; then
    echo "Static input directory not found: $static_dir" >&2
    exit 1
  fi

  mkdir -p "$output_dir"

  # Normalize base paths to avoid duplicate slashes later.
  scratch_base="${scratch_base%/}"
  repo_root="${repo_root%/}"

  local generated=0
  while IFS= read -r -d '' input_file; do
    local filename stem rel_path remote_input job_file
    filename="$(basename "$input_file")"
    # Skip base templates or non-run input decks explicitly.
    if [[ "$filename" == *"_basefile"*.i ]]; then
      continue
    fi

    stem="${filename%.i}"
    rel_path="${input_file#"$repo_root/"}"
    remote_input="$scratch_base/$rel_path"
    job_file="$output_dir/job_frontera_${stem}.sbatch"

    cat > "$job_file" <<SBATCH
#!/bin/bash
#SBATCH -J ${stem}        # Job name
#SBATCH -o ${stem}.o%j    # Name of stdout output file
#SBATCH -e ${stem}.e%j    # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 15              # Total # of nodes 
#SBATCH -n 400             # Total # of mpi tasks
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A EAR20006        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=chunhui3@illinois.edu

# Load necessary modules
module swap intel gcc
#module swap impi mvapich2-x
module load cuda
export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Set compilers and flags
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77

export CXXFLAGS=-I/opt/apps/gcc/9.1.0/include/c++/9.1.0/
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77

# Enable MPI debugging
export MV2_DEBUG=1
export MV2_SHOW_ENV_INFO=1

export MOOSE_JOBS=6 METHODS=opt

ibrun ${executable} -i ${remote_input} --allow-unused
SBATCH

    chmod +x "$job_file"
    echo "Wrote $job_file"
    ((generated++))
  done < <(find "$static_dir" -maxdepth 1 -type f -name 'static_solve_*.i' -print0)

  if [[ "$generated" -eq 0 ]]; then
    echo "No sbatch files generated (check for input decks or applied filters)."
  else
    echo "Done. Generated $generated sbatch files in $output_dir."
  fi
}

main "$@"
