#!/bin/bash
set -e

REPO="$(git rev-parse --show-toplevel)"
MIRROR="$(dirname "$REPO")/PULSUS_Manuscript_Overleaf"
MSG="${1:-Sync manuscript from main repository}"

echo "Updating Overleaf manuscript worktree..."
git -C "$MIRROR" pull --ff-only

echo
echo "Syncing paper/ -> Overleaf manuscript worktree..."

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
  "$REPO/paper/" \
  "$MIRROR/"

echo

if git -C "$MIRROR" diff --quiet && \
   git -C "$MIRROR" diff --cached --quiet && \
   [ -z "$(git -C "$MIRROR" ls-files --others --exclude-standard)" ]; then
    echo "No manuscript changes to push."
    exit 0
fi

echo "Committing manuscript changes..."
git -C "$MIRROR" add -A
git -C "$MIRROR" commit -m "$MSG"

echo
echo "Pushing manuscript to Overleaf..."
git -C "$MIRROR" push overleaf-paper HEAD:main

echo
echo "Manuscript synchronized with Overleaf."
