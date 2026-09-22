"""blaKPC allele typing and ceftazidime-avibactam marker interpretation.

Ceftazidime-avibactam resistance in KPC-producing Enterobacterales is driven
overwhelmingly by amino-acid changes in the KPC enzyme itself.  Two independent
lines of evidence are used here, in this order:

1. **Allele identity.**  The observed set of differences from KPC-2 is matched
   against every blaKPC allele defined by NCBI (``data/blaKPC_alleles.tsv``).
   An exact match names the allele and carries NCBI's own curation of whether
   that allele is an inhibitor-resistant, extended-spectrum enzyme.

2. **Mechanism.**  Changes are checked against the structural hotspots that the
   literature ties to avibactam escape - the Omega loop and the 237-243 region
   around the active site, plus in-frame insertions/duplications in the
   266-276 loop.

A change that is in a hotspot but is not a documented variant yields an
*Indeterminate* call rather than a Resistant one: the tool says it does not
know, instead of guessing.

References
----------
* Ambler numbering and the absent positions 58/253: Standardized numbering and
  alignment of the KPC family of beta-lactamases, Antimicrob Agents Chemother
  (2025), doi:10.1128/aac.01868-25
* D179 substitutions and Omega-loop destabilisation: Antimicrob Agents
  Chemother (2022), doi:10.1128/aac.02414-21
* Variant overview, incl. V240G, T243M, Omega-loop indels and the 266-276
  insertions: Klebsiella pneumoniae Carbapenemase Variants Resistant to
  Ceftazidime-Avibactam: an Evolutionary Overview, Antimicrob Agents Chemother
  (2022), doi:10.1128/aac.00447-22
"""

import csv
import re
from pathlib import Path

_DATA_DIR = Path(__file__).parent / 'data'
_ALLELE_TABLE = _DATA_DIR / 'blaKPC_alleles.tsv'

# Hotspot regions, in Ambler numbering.
OMEGA_LOOP = range(164, 180)          # 164-179
ACTIVE_SITE_237_243 = range(237, 244)  # 237-243
INSERTION_LOOP_266_276 = range(266, 277)

# Substitutions with direct, repeatedly published evidence of conferring
# ceftazidime-avibactam resistance.
DOCUMENTED_SUBSTITUTIONS = {
    'D179Y', 'D179N', 'D179G', 'D179V',
    'L169P', 'L169Q', 'L169M',
    'V240G', 'V240A',
    'T243M', 'T243A',
    'E166K',
}

_CHANGE_PATTERN = re.compile(r'^([A-Z])(\d+)(?:_([A-Z])(\d+))?(del|[A-Z])$')
_INSERTION_PATTERN = re.compile(r'^ins([A-Z]+)@(\d+)$')


def _load_allele_table():
    alleles = {}
    if not _ALLELE_TABLE.exists():
        return alleles
    with open(_ALLELE_TABLE) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = _normalise_change_set(row['Changes_Ambler'].split(';'))
            alleles[key] = row
    return alleles


def _normalise_change_set(labels):
    return tuple(sorted(label for label in labels if label))


_ALLELES = _load_allele_table()


def identify_kpc_allele(change_labels):
    """Name the blaKPC allele whose definition is exactly this change set.

    Returns the allele table row, or ``None`` when the combination is not a
    described allele (a novel variant).
    """
    return _ALLELES.get(_normalise_change_set(change_labels))


def _positions_touched(label):
    """Every Ambler position a change label refers to."""
    insertion = _INSERTION_PATTERN.match(label)
    if insertion:
        return [int(insertion.group(2))], 'insertion'

    match = _CHANGE_PATTERN.match(label)
    if not match:
        return [], 'unknown'

    start = int(match.group(2))
    end = int(match.group(4)) if match.group(4) else start
    kind = 'deletion' if match.group(5) == 'del' else 'substitution'
    return list(range(start, end + 1)), kind


def assess_kpc_changes(change_labels, allele_row=None):
    """Interpret KPC changes for ceftazidime-avibactam.

    Returns ``{'call': 'Resistant'|'Indeterminate'|'Susceptible', 'evidence': [...]}``
    where the call describes the contribution of this enzyme only.
    """
    evidence = []
    confident = False
    uncertain = False

    if allele_row and allele_row.get('Inhibitor_resistant') == 'yes':
        confident = True
        evidence.append(
            f"{allele_row['Allele']} is curated by NCBI as "
            f"\"{allele_row['NCBI_product_name']}\" "
            f"(subclass {allele_row['NCBI_subclass'] or 'unassigned'})"
        )

    for label in change_labels:
        if not label:
            continue
        positions, kind = _positions_touched(label)

        if label in DOCUMENTED_SUBSTITUTIONS:
            confident = True
            evidence.append(f"{label}: documented ceftazidime-avibactam resistance substitution")
            continue

        in_omega = any(position in OMEGA_LOOP for position in positions)
        in_active_site = any(position in ACTIVE_SITE_237_243 for position in positions)
        in_insertion_loop = any(position in INSERTION_LOOP_266_276 for position in positions)

        if kind in ('deletion', 'insertion') and in_omega:
            confident = True
            evidence.append(
                f"{label}: in-frame {kind} in the Omega loop (Ambler 164-179), "
                f"the best-described route to avibactam escape")
        elif in_omega:
            uncertain = True
            evidence.append(
                f"{label}: substitution in the Omega loop (Ambler 164-179) that is not "
                f"a documented variant - effect unknown")
        elif in_active_site:
            uncertain = True
            evidence.append(
                f"{label}: change in the 237-243 active-site region that is not a "
                f"documented variant - effect unknown")
        elif kind == 'insertion' and in_insertion_loop:
            uncertain = True
            evidence.append(
                f"{label}: in-frame insertion in the 266-276 loop, a region where "
                f"insertions have been reported in ceftazidime-avibactam-resistant alleles")

    if confident:
        return {'call': 'Resistant', 'evidence': evidence}
    if uncertain:
        return {'call': 'Indeterminate', 'evidence': evidence}
    return {'call': 'Susceptible', 'evidence': evidence}


def describe_kpc_result(change_labels, allele_row):
    """One-line description of what the KPC copy is."""
    if allele_row:
        return allele_row['Allele']
    if not change_labels:
        return 'blaKPC-2'
    return 'novel blaKPC variant (' + ', '.join(change_labels) + ')'
