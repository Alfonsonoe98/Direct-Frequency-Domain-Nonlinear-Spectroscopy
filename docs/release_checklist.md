# Preparing the plotting release

The plotting additions are currently a development candidate, version `0.3.0.dev0`. They have not been published to PyPI.

Before publishing 0.3.0:

1. Review the plotting code, tests, examples, and manual together.
2. Change the development version to `0.3.0` in package metadata and documentation, remove the development-install notices, and rebuild the reviewed PDF.
3. Build distributions, install the wheel in a clean environment, and run the tests and documented examples against that installed wheel.
4. Confirm the PyPI Trusted Publisher uses the renamed GitHub repository, `release.yml`, and the `pypi` environment.
5. Commit and push the reviewed release changes, then publish a GitHub release tagged `v0.3.0`. The existing workflow publishes the tagged package to PyPI.
6. Verify the PyPI installation and Zenodo archive; update citation metadata with the actual new archive identifier. Do not relabel the archived 0.2.0 DOI as the 0.3.0 release.

`CITATION.cff` currently refers to the archived 0.2.0 release. Keep that association until a new archive identifier is available.
