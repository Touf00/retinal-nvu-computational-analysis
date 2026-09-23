# GitHub-to-Zenodo release checklist

Do not create the archival release until the cleaned notebook has completed from raw public inputs and every assertion has passed.

## Before the first release

1. Replace the entity placeholder in `CITATION.cff` with the real software authors, ORCIDs, and affiliations. Validate the file with the [Citation File Format validator](https://citation-file-format.github.io/cff-initializer-javascript/#/validation).
2. Confirm that every contributor accepts the MIT license.
3. Run `python tests/test_release_package.py`.
4. Run the notebook top to bottom in the declared environment. Confirm that it ends with `All fail-fast validation checks passed.`
5. Review `results/generated/logs/input_manifest_sha256.csv`, the generated summary tables, and both figure formats.
6. Commit the verified outputs you actually want archived. Keep raw inputs and large cell-level intermediates excluded.

## Connect the GitHub repository to Zenodo

1. Sign in to [Zenodo](https://zenodo.org/) and link the correct GitHub account if it is not already connected.
2. Open the profile menu, choose **GitHub**, and click **Sync now**.
3. Find `retinal-nvu-computational-analysis` and enable its repository toggle.
4. Refresh the Zenodo GitHub page and verify that the repository is enabled.

## Mint the first version DOI

1. On GitHub, create a release from the verified commit. Use a semantic tag such as `v1.0.0` and a clear release title.
2. Publish the GitHub release. Zenodo should ingest the release automatically after the repository is enabled.
3. Wait for processing, then open the record DOI from the Zenodo GitHub integration page.
4. Check the record title, creators, description, version, license, keywords, related identifiers, and uploaded archive before treating the DOI as final metadata.
5. Under **External resources**, verify the archival status. Zenodo may take additional time to complete downstream archival.
6. Add the version DOI badge to `README.md`. For manuscripts and software citation, cite the version DOI; use the concept DOI when referring to the evolving software project across releases.
7. Add the repository URL and DOI to `CITATION.cff`, commit that metadata update, and create a later patch release if you need the revised citation file itself archived.

Zenodo’s current official instructions are at:

- https://help.zenodo.org/docs/github/enable-repository/
- https://help.zenodo.org/docs/github/archive-software/github-upload/

## Failure modes worth checking

- Zenodo cannot see the repository: unlink/relink the GitHub integration, click **Sync now**, and check organization access.
- Release ingestion fails: open the failed release entry in Zenodo and inspect **Errors**, especially malformed citation metadata.
- Wrong authors appear: fix `CITATION.cff` before publishing another GitHub release; do not paper over incorrect authorship in prose.
- Huge archive or raw data appear in the release: stop and remove them from Git history before publishing the archival tag.
