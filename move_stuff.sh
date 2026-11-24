#!/bin/bash
set -euo pipefail

BASE="/projects/b1169/boles/als_cns_visium"
cd "$BASE"

mkdir -p data/07_c2l

COMBOS=("mcx_discrete" "mcx_microglia_subtypes")

for combo in "${COMBOS[@]}"; do
  old="results/07_c2l/${combo}"

  echo "=== $combo ==="

  # Data objects -> data/07_c2l/, flat and combo-prefixed
  mv -v "${old}/models/model" "data/07_c2l/${combo}_regression_model"
  mv -v "${old}/adata_objects/ref_after_training.h5ad" "data/07_c2l/${combo}_ref_after_training.h5ad"

  # Diagnostic plots/CSVs -> flatten up out of plots/ into the combo results folder
  mv -v "${old}"/plots/* "${old}/"

  # Drop the old inf_aver.csv -- 07a no longer writes this (07b recomputes an
  # identical table and writes its own copy once it runs), so this is the
  # redundant file dropped in the directory rework, not something the new
  # script layout expects to see here.
  rm -fv "${old}/inf_aver.csv"

  # Now-empty old subfolders
  rmdir "${old}/plots" "${old}/models" "${old}/adata_objects"
done

echo "Done. New layout:"
find data/07_c2l results/07_c2l -maxdepth 2 | sort