#!/usr/bin/env bash
# Run from the app folder, AFTER the current analysis finishes and Streamlit stops.
set -euo pipefail
if [ ! -f app.py ] || [ ! -d nanopore ]; then
  echo 'Open the terminal in your app folder (the folder containing app.py) first.' >&2
  exit 1
fi
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
base='https://raw.githubusercontent.com/monikasafhauzer/project/nanopore-race-app'
files=(app.py requirements.txt nanopore/alignment.py nanopore/summary.py nanopore/workflow.py nanopore/jobs.py nanopore/worker.py nanopore/job_ui.py nanopore/dal2.py nanopore/viewer.py nanopore/viewer.html pages/1_DAL2_Alignment_Viewer.py)
for file in "${files[@]}"; do
  mkdir -p "$stage/$(dirname "$file")"
  curl -fL "$base/$file" -o "$stage/$file"
done
python_bin="$HOME/.nanopore-app/bin/python"
uv_bin="$HOME/.nanopore-tools/bin/uv"
if [ ! -x "$python_bin" ] || [ ! -x "$uv_bin" ]; then
  echo 'Your prepared ~/.nanopore-app and ~/.nanopore-tools environments were not found. See README.' >&2
  exit 1
fi
"$uv_bin" pip install --python "$python_bin" -r "$stage/requirements.txt"
backup="$(pwd)/update-backups/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$backup"
for file in "${files[@]}"; do
  if [ -f "$file" ]; then
    mkdir -p "$backup/$(dirname "$file")"
    cp "$file" "$backup/$file"
  fi
  mkdir -p "$(dirname "$file")"
  cp "$stage/$file" "$file"
done
echo "Updated. Previous files are backed up in $backup"
echo 'Start the app with: ~/.nanopore-app/bin/python -m streamlit run app.py'
