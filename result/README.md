# Recorded application

These files transcribe the author's reported MCH, MCV and RDW production outcomes and prefix-pilot results. The author reported successful final integrity auditing on 14 November 2026. They are not reruns of the full datasets in the release-building environment.

- `production_qc.tsv`: original production summary columns and reported values.
- `production_qc.json`: the same summary plus counters for automated conservation tests. PLUS_PASS and MINUS_PASS are uniquely reconstructed from the reported integer totals and rounded strand-specific rate; that derivation is recorded in the JSON.
- `pilot_qc.tsv`: first 10,000 rows per trait; all passed on the positive strand. This prefix sample is not a genome-wide success-rate estimate.

Across traits, 150,741,969 input records yielded 141,659,492 accepted SNP records and 9,082,477 review records. The 9,001,714 records marked non-SNP/unsupported allele remain outside the accepted SNP set. Record counts overlap across traits and do not represent that many distinct loci.

The actual full input data, run directories, binary fingerprints and absolute-path manifests are retained with the author's completed run. They are not reconstructed or invented here. The provenance exporter can be used on that run to prepare a path-redacted metadata supplement.
