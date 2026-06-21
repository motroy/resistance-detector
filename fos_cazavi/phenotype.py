"""
Genotype-to-phenotype prediction for fosfomycin (FOS) and
ceftazidime-avibactam (CAZ/AVI).

These calls are derived purely from the genes/mutations detected by the
BLAST and GAMMA/SeqKit pipelines, using literature-established resistance
mechanisms (see README/docs for the supporting references). They are a
genotypic prediction, NOT a substitute for phenotypic antimicrobial
susceptibility testing (AST).
"""

# Acquired, plasmid-borne fosA-family enzymes: presence alone inactivates
# fosfomycin and is sufficient for a resistant call.
FOS_ACQUIRED_GENES = {'fosA', 'fosA3', 'fosA4', 'fosA5', 'fosA7', 'fosA11'}

# Chromosomal fosfomycin uptake/regulatory genes: a loss-of-function (premature
# stop codon) mutation reduces fosfomycin uptake and confers resistance. The
# near-universal fosAKP I91V variant is intentionally excluded - it is not a
# loss-of-function change and is present in genotypically-susceptible strains.
# unified_results normalizes gene names to lowercase (see
# MutationDetector._normalize_gene_name in fos_cazavi/mutations.py), so this
# set is matched case-insensitively.
FOS_TRANSPORT_GENES = {
    g.lower() for g in
    ('murA', 'uhpT', 'uhpA', 'uhpB', 'uhpC', 'glpT', 'cyaA', 'ptsI', 'galU', 'lon')
}

# Tracked blaKPC Omega-loop/X-loop substitutions known to confer
# ceftazidime-avibactam resistance (see fos_cazavi/utils.py KNOWN_MUTATIONS).
KPC_AVIBACTAM_RESISTANCE_MARKERS = (
    'L166', 'E167', 'L168', 'N169', 'D179', 'V240', 'T243'
)


def _is_lof_mutation(mutation_str):
    """A mutation name/token represents a loss-of-function change if it
    encodes a premature stop codon ('*') or an in-frame deletion ('del')."""
    return '*' in mutation_str or 'del' in mutation_str


def _is_kpc_avibactam_marker(mutation_str):
    """A mutation name/token represents a tracked Omega-loop/X-loop
    substitution or indel known to confer avibactam resistance."""
    if 'del' in mutation_str or 'ins' in mutation_str:
        return True
    return any(marker in mutation_str for marker in KPC_AVIBACTAM_RESISTANCE_MARKERS)


def predict_fos_phenotype(blast_results, unified_results):
    """Predict fosfomycin susceptibility from detected genes/mutations.

    Returns a dict: {'phenotype': 'Resistant'|'Susceptible', 'evidence': [str, ...]}
    """
    evidence = []

    for r in blast_results or []:
        if r['gene'] in FOS_ACQUIRED_GENES:
            evidence.append(
                f"Acquired fosfomycin-inactivating enzyme {r['gene']} detected "
                f"({r['identity']}% identity, {r['coverage']}% coverage)"
            )

    for r in unified_results or []:
        if r['gene'].lower() in FOS_TRANSPORT_GENES and _is_lof_mutation(r['mutation']):
            evidence.append(
                f"Loss-of-function mutation {r['mutation']} in fosfomycin "
                f"uptake/regulatory gene {r['gene']} (reduced fosfomycin uptake)"
            )

    phenotype = 'Resistant' if evidence else 'Susceptible'
    if not evidence:
        evidence.append(
            "No acquired fosA-family enzyme or loss-of-function mutation in "
            "fosfomycin uptake/regulatory genes detected"
        )

    return {'phenotype': phenotype, 'evidence': evidence}


def predict_cazavi_phenotype(blast_results, unified_results):
    """Predict ceftazidime-avibactam susceptibility from detected genes/mutations.

    Returns a dict: {'phenotype': 'Resistant'|'Susceptible', 'evidence': [str, ...]}
    """
    evidence = []

    for r in blast_results or []:
        if 'KPC' not in r['gene'].upper():
            continue
        if r['mutations'] != '-' and _is_kpc_avibactam_marker(r['mutations']):
            evidence.append(
                f"blaKPC variant {r['gene']} carries Omega-loop/X-loop "
                f"resistance marker(s): {r['mutations']}"
            )

    for r in unified_results or []:
        if 'KPC' not in r['gene'].upper():
            continue
        if _is_kpc_avibactam_marker(r['mutation']):
            evidence.append(
                f"blaKPC mutation {r['mutation']} detected in {r['gene']} "
                f"(Omega-loop/X-loop avibactam-resistance marker)"
            )

    phenotype = 'Resistant' if evidence else 'Susceptible'
    if not evidence:
        evidence.append(
            "No blaKPC Omega-loop/X-loop resistance marker detected "
            "(wildtype blaKPC, if present, remains avibactam-inhibitable)"
        )

    return {'phenotype': phenotype, 'evidence': evidence}


def predict_phenotypes(blast_results, unified_results):
    """Predict both FOS and CAZ/AVI phenotypes.

    Returns:
        {
          'fosfomycin': {'phenotype': ..., 'evidence': [...]},
          'ceftazidime_avibactam': {'phenotype': ..., 'evidence': [...]},
          'disclaimer': str
        }
    """
    return {
        'fosfomycin': predict_fos_phenotype(blast_results, unified_results),
        'ceftazidime_avibactam': predict_cazavi_phenotype(blast_results, unified_results),
        'disclaimer': (
            "Genotype-based prediction only; not a substitute for phenotypic "
            "antimicrobial susceptibility testing (AST)."
        )
    }
