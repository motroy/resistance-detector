# PRJNA741867 — validation against a published clinical study

Six *Klebsiella pneumoniae* ST307 assemblies from BioProject PRJNA741867. The
study reported three paired isolates: a ceftazidime-avibactam-susceptible
baseline and a resistant isolate that emerged on therapy, in each of three
patients (A, B, C).

The pipeline was run blind to those labels:

```bash
fos-cazavi fos-cazavi-all -a <assembly>.fna \
    -d fos_cazavi/data/example_database.fasta \
    -o <sample> --organism Klebsiella_pneumoniae
```

## Results

| Sample | Accession | Paper: CAZ/AVI | blaKPC allele called | Changes (Ambler) | Predicted CAZ/AVI | Agreement |
|---|---|---|---|---|---|---|
| A-1 (S16)   | GCA_022423605.1 | Susceptible | blaKPC-3  | H274Y | Susceptible | ✅ |
| A-2 (1-A3)  | GCA_022423565.1 | Resistant   | blaKPC-46 | L169P; H274Y | Resistant | ✅ |
| B-1 (1-G9)  | GCA_022423575.1 | Susceptible | blaKPC-3  | H274Y | Susceptible | ✅ |
| B-2 (2-C6)  | GCA_022423525.1 | Resistant   | blaKPC-66 | E166_L167del; H274Y | Resistant | ✅ |
| C-1 (1-I6)  | GCA_022423535.1 | Susceptible | blaKPC-3  | H274Y | Susceptible | ✅ |
| C-2 (HC)    | GCA_022423505.1 | Resistant   | blaKPC-92 | E168D; L169_N170del; H274Y | Resistant | ✅ |

6/6 concordant, including exact allele assignment for all three resistant
isolates.

## Reading the numbering

Positions are given in standardised Ambler numbering, which omits positions 58
and 253. Papers that number KPC sequentially without allowing for those absent
positions report the same variants one lower — what is written here as
`L169P` (KPC-46) appears as `L168P` in some publications, and `E166_L167del`
(KPC-66) as `E167_L168del`. The sequence change is identical; only the label
differs. See doi:10.1128/aac.01868-25.

`H274Y` is carried by all six isolates: it is what distinguishes KPC-3 from
KPC-2 and is **not** an avibactam-escape mutation. The pipeline lists it among
the protein changes but does not treat it as a resistance marker, which is why
the three baseline isolates are called Susceptible.

## Fosfomycin

All six were predicted fosfomycin-susceptible. Each carries the intrinsic
chromosomal `fosAKP`, which is present in fosfomycin-susceptible *K.
pneumoniae* and is deliberately not scored as acquired resistance; no acquired
fosA-family enzyme and no loss-of-function change in the uptake/regulatory
genes was found.

## What this does and does not show

It shows that the blaKPC typing and the ceftazidime-avibactam logic reproduce
the published phenotypes on real clinical assemblies, including two indel
variants. It does not validate the fosfomycin side (no resistant isolate here),
and six genomes from one clonal background is a narrow test — see
[../../docs/METHODS.md](../../docs/METHODS.md) for the method's limits.
