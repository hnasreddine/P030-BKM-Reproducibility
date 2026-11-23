# GitHub public + Zenodo version v1.0.0

1. Confirm that this is a **new public repository**, next to (not inside)
   `hnasreddine/P014-Tetrahedral-Quantum-Certificate`, and contains neither
   the unpublished article nor AI/internal audit reports.
2. Install the tested environment: `python -m pip install -r requirements-tested.txt`.
3. Run `python tools/verify_release.py`, then `bash RUN_PRECOMPUTED.sh`.
4. Run `bash RUN_EXACT.sh` on a supported IEEE-754 Linux x86-64 machine;
   review each mathematical source in `docs/PROOF_MAP.md`.
5. Commit all release files to the separate public GitHub repository;
   inspect GitHub Actions run on that exact commit. Do not publish a tag or
   release while CI is red, or before reviewing rights and licensing.
6. Enable the repository in Zenodo GitHub integration, **before** creating
   the immutable `v1.0.0` GitHub release. Tag the exact passing commit.
7. Verify the Zenodo deposit content, assigned DOI and GitHub release link.
   Add that real version DOI to `CITATION.cff` only in a subsequent, carefully
   versioned change; never mutate the archived `v1.0.0` snapshot.
8. Update the private P030 JMP manuscript *Data Availability* statement
   with the actual DOI and access date, not a placeholder or invented identifier.

The current package is scientifically packaged and content-verified; a public
GitHub repository and DOI are NOT claimed to exist until verified externally.
