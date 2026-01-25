# Test Notes

## Database Creation Issue
The `create_reference_database.py` script downloads full GenBank records for `fosA4` (KP324830.1) and `fosA5` (KY270852.1). These accessions correspond to complete plasmids (82kb and 282kb respectively), not just the resistance genes.

The `resistance_detector.py` calculates coverage as:
```python
coverage = (length / slen) * 100
```
where `slen` is the subject length (length of the sequence in the database).

Because the database contains the full plasmids, the subject length is very large. A perfect match of the gene (approx 500bp) results in < 1% coverage.
`500 / 82000 = 0.6%`

Since the default minimum coverage is 80%, `fosA4` and `fosA5` will **never be detected** using the database created by the standard script.

## Test Results
- **E. coli K-12 MG1655 (NC_000913.3)**: No resistance genes detected. This is expected as it is a susceptible lab strain.
- **K. pneumoniae NTUH-K2044 (GCA_000009885.1)**: No resistance genes detected. This strain contains chromosomal *fosA* (homologous to *fosA* variants), but it was not detected. This could be due to:
    1. Divergence from the specific *fosA3/4/5* variants in the database.
    2. The database issue described above (if it matched `fosA5` plasmid backbone or gene but failed coverage).

## Recommendations
The `create_reference_database.py` script should be updated to extract only the coding sequences for the resistance genes, or the accessions should be changed to point to gene-specific entries if available.
