#!/usr/bin/env bash
# Show scheduler state together with committed data results and recent log output.
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
config_path="${1:-$repo_root/config.example.json}"
python_bin="${GWAMA_PYTHON:-$(command -v python3)}"
if [[ -f "$repo_root/logs/last_job_id.txt" ]]; then
  read -r job_id < "$repo_root/logs/last_job_id.txt"
  if [[ "$job_id" =~ ^[0-9]+$ ]]; then
    if command -v squeue >/dev/null; then
      squeue -j "$job_id" -o "%.18i %.12P %.20j %.10T %.6C %.12M %.30R"
    fi
    if command -v sacct >/dev/null; then
      sacct -j "$job_id" --format=JobID,Partition,State,Elapsed,AllocCPUS,MaxRSS
    fi
  fi
fi
"$python_bin" -u "$repo_root/scripts/Step_03_liftover.py" status --config "$config_path"
