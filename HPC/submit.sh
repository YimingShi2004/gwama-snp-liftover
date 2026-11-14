#!/usr/bin/env bash
# Submit one sequential multi-trait production job.
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
config_path="${1:-$repo_root/config.example.json}"
if [[ "$config_path" != /* ]]; then config_path="$PWD/$config_path"; fi
python_bin="${GWAMA_PYTHON:-$(command -v python3)}"
partition="${GWAMA_PARTITION:-bi2}"
mkdir -p "$repo_root/logs"
printf -v wrapped 'env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 %q -u %q run --config %q' \
  "$python_bin" "$repo_root/scripts/Step_03_liftover.py" "$config_path"
if job_id="$(sbatch --parsable --partition="$partition" --job-name=GWAMA_hg38 \
  --cpus-per-task=4 --mem=24G --time=24:00:00 --chdir="$repo_root" \
  --output="$repo_root/logs/Step_03_GWAS_hg38_%j.log" --wrap="$wrapped")"; then
  printf '%s\n' "${job_id%%;*}" > "$repo_root/logs/last_job_id.txt"
  printf 'Submitted job %s to %s\n' "$job_id" "$partition"
fi
