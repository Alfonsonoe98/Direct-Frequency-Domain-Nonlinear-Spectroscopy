#!/bin/bash
set -e

REPO="$(git rev-parse --show-toplevel)"

cd "$REPO"

if [ -n "$(git status --porcelain)" ]; then
    echo "Working tree is not clean. Commit or stash your changes first."
    exit 1
fi

echo "Fetching latest Technical Notes Overleaf history..."
git fetch overleaf-notes

echo
echo "Checking for Overleaf changes..."
git merge --ff-only overleaf-notes/main

echo
echo "Pushing repository to Technical Notes Overleaf..."
git push overleaf-notes HEAD:main

echo
echo "Technical Notes Overleaf synchronized."
