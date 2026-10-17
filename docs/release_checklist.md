# Publishing PULSUS 0.3.0

The plotting release is prepared in [pull request #1](https://github.com/Alfonsonoe98/PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators/pull/1), from `feature/linear-plotting` into `main`. Package metadata must say `0.3.0` before tagging. A merge into `main` does **not** publish to PyPI: publishing the GitHub release triggers `.github/workflows/release.yml`.

## Before merging

1. Review the plotting code, tests, standalone examples, shortened README, API reference, and rebuilt PDF together.
2. Check that the complete test suite passes, distributions build, and the installed wheel runs the documented examples. Record the actual checks in [VALIDATION.md](VALIDATION.md); scientific benchmarks are described separately in [method_notes.md](method_notes.md).
3. Check that PyPI's Trusted Publisher has these exact settings:

   | Setting | Value |
   |---|---|
   | GitHub owner | `Alfonsonoe98` |
   | Repository | `PULSUS-Pulse-dressed-Liouvillian-Spectroscopy-Using-Superoperators` |
   | Workflow filename | `release.yml` |
   | Environment | `pypi` |

4. Commit and push the reviewed preparation changes to `feature/linear-plotting`. Mark PR #1 ready for review if it is still a draft (`gh pr ready 1`).

## Merge into main

Run these commands from the local repository. The branch commands preserve the untracked `prototypes/` directory; they do not clean or remove untracked files.

```bash
cd "/Users/alfonsocastillogonzalez/Desktop/Research/Spectroscopy_proyect/GitHub/PULSUS"
gh pr checks 1 --watch
git switch main
gh pr merge 1 --merge --delete-branch
git pull --ff-only origin main
```

This merges the reviewed PR, deletes its feature branch, and brings local `main` to the merged commit. Inspect `git status --short` and `git log -1 --oneline`; an untracked `prototypes/` directory may remain. Do not include it in the release unless deliberately reviewed and added.

## Publish the release

After merging and updating local `main`, create the version tag at that commit and publish the reviewed release notes:

```bash
git tag -a v0.3.0 -m "PULSUS v0.3.0"
git push origin v0.3.0
gh release create v0.3.0 --verify-tag --title "PULSUS v0.3.0" --notes-file docs/releases/0.3.0.md
```

Publishing the GitHub release starts **Publish to PyPI** in the repository's Actions tab. Wait for both the build and publish jobs to succeed before announcing availability. Do not use the workflow's manual dispatch to publish the same version again.

## After publication

1. Check the [PyPI release page](https://pypi.org/project/pulsus-spec/0.3.0/) and install the published package in a fresh environment with `python -m pip install "pulsus-spec==0.3.0"`.
2. Run a linear and a two-dimensional plotting example from that installed package, including optional figure export.
3. Check the actual Zenodo archive for the new release before updating citation metadata. `CITATION.cff` currently identifies the archived 0.2.0 release; keep its DOI associated with 0.2.0 until the new archive identifier is available.

The tag and release refer to the reviewed `main` commit. Subsequent documentation or citation corrections should be committed normally without moving that tag.
