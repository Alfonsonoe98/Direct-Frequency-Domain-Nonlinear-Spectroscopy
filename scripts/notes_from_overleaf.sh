#!/bin/bash
set -e

REPO="$(git rev-parse --show-toplevel)"

cd "$REPO"

if [ -n "$(git status --porcelain)" ]; then
    echo "Working tree is not clean. Commit or stash your changes first."
    exit 1
fi

echo "Fetching latest Technical Notes changes from Overleaf..."
git fetch overleaf-notes

echo
echo "Updating local branch..."
git merge --ff-only overleaf-notes/main

echo
echo "Technical Notes changes synchronized locally."
echo "Remember to push them to GitHub with: git push"
