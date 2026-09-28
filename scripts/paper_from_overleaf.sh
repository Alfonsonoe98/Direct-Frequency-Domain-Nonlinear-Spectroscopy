#!/bin/bash
set -e

REPO="$(git rev-parse --show-toplevel)"
MIRROR="$(dirname "$REPO")/PULSUS_Manuscript_Overleaf"

echo "Pulling latest manuscript from Overleaf..."
git -C "$MIRROR" pull --ff-only

echo
echo "Syncing Overleaf manuscript -> paper/..."

rsync -av --delete \
  --exclude='.git' \
  --exclude='*.aux' \
  --exclude='*.bbl' \
  --exclude='*.blg' \
  --exclude='*.bcf' \
  --exclude='*.fdb_latexmk' \
  --exclude='*.fls' \
  --exclude='*.log' \
  --exclude='*.out' \
  --exclude='*.run.xml' \
  --exclude='*.synctex.gz' \
  --exclude='main.pdf' \
  "$MIRROR/" \
  "$REPO/paper/"

echo
echo "Changes in main repository paper/:"
git -C "$REPO" status --short -- paper
