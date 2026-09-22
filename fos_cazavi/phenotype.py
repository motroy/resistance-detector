"""Genotype-to-phenotype prediction for fosfomycin and ceftazidime-avibactam.

Every call is one of three values:

``Resistant``
    A mechanism with established, published evidence is present.
``Indeterminate``
    Something relevant was found but its effect is not established - a novel
    change in a resistance hotspot, or a target gene that could not be assessed
    because it ran off a contig.  The tool says so instead of defaulting to
    susceptible.
``Susceptible``
    The known mechanisms were looked for and not found.

The result is a genotypic prediction and never a substitute for phenotypic
antimicrobial susceptibility testing.
"""

from . import betalactamase
from .references import (
    FOSA_FAMILIES, FOS_TRANSPORT_GENES, INTRINSIC_GENES, gene_family, is_mbl,
)

RESISTANT = 'Resistant'
INDETERMINATE = 'Indeterminate'
SUSCEPTIBLE = 'Susceptible'

DISCLAIMER = ('Genotype-based prediction only; not a substitute for phenotypic '
              'antimicrobial susceptibility testing (AST).')


def _resolve(resistant_evidence, uncertain_evidence, nothing_found_message):
    if resistant_evidence:
        return {'phenotype': RESISTANT,
                'evidence': resistant_evidence + uncertain_evidence}
    if uncertain_evidence:
        return {'phenotype': INDETERMINATE, 'evidence': uncertain_evidence}
    return {'phenotype': SUSCEPTIBLE, 'evidence': [nothing_found_message]}


def predict_fos_phenotype(blast_results, unified_results=None):
    """Predict fosfomycin susceptibility."""
    resistant, uncertain = [], []

    for result in blast_results or []:
        gene = result['gene']
        family = gene_family(gene)

        if gene in INTRINSIC_GENES:
            # Intrinsic chromosomal fosA (e.g. fosAKP in K. pneumoniae) is
            # present in susceptible isolates and is not scored here.
            continue

        if family in FOSA_FAMILIES:
            resistant.append(
                f"Acquired fosfomycin-modifying enzyme {result['allele']} "
                f"({result['identity']}% identity, {result['coverage']}% coverage)")
            continue

        if gene in FOS_TRANSPORT_GENES or gene == 'murA':
            if result.get('loss_of_function'):
                if result.get('complete'):
                    resistant.append(
                        f"Loss of function in fosfomycin uptake/regulatory gene {gene}: "
                        f"{result['lof_description']}")
                else:
                    uncertain.append(
                        f"Apparent loss of function in {gene} "
                        f"({result['lof_description']}), but the gene runs off a contig "
                        f"boundary - could be an assembly artefact")
            elif result.get('reported_mutations'):
                resistant.append(
                    f"Curated fosfomycin-resistance mutation(s) in {gene}: "
                    f"{', '.join(result['reported_mutations'])}")

    return _resolve(
        resistant, uncertain,
        'No acquired fosfomycin-modifying enzyme and no loss-of-function change '
        'in the fosfomycin uptake/regulatory genes were detected')


def predict_cazavi_phenotype(blast_results, unified_results=None):
    """Predict ceftazidime-avibactam susceptibility."""
    resistant, uncertain = [], []
    kpc_seen = False

    for result in blast_results or []:
        gene = result['gene']
        family = gene_family(gene)

        if is_mbl(gene):
            resistant.append(
                f"Metallo-beta-lactamase {result['allele']} detected "
                f"({result['identity']}% identity) - avibactam does not inhibit "
                f"metallo-enzymes, so ceftazidime-avibactam is not active")
            continue

        if family != 'blaKPC':
            continue

        kpc_seen = True
        assessment = betalactamase.assess_kpc_changes(
            result.get('changes', []), result.get('allele_row'))
        prefix = f"{result['allele']} on {result['contig']}"

        if not result.get('complete'):
            uncertain.append(
                f"{prefix}: gene runs off a contig boundary, so its variant "
                f"content could not be fully assessed")

        if assessment['call'] == 'Resistant':
            resistant.extend(f"{prefix}: {item}" for item in assessment['evidence'])
        elif assessment['call'] == 'Indeterminate':
            uncertain.extend(f"{prefix}: {item}" for item in assessment['evidence'])

    if kpc_seen and not resistant and not uncertain:
        nothing_found = ('blaKPC detected but carrying no Omega-loop, 237-243 or '
                         'insertion-loop change associated with avibactam escape; '
                         'avibactam is expected to inhibit it')
    else:
        nothing_found = ('No ceftazidime-avibactam resistance mechanism detected '
                         '(no metallo-beta-lactamase, no blaKPC escape variant)')

    return _resolve(resistant, uncertain, nothing_found)


def predict_phenotypes(blast_results, unified_results=None):
    return {
        'fosfomycin': predict_fos_phenotype(blast_results, unified_results),
        'ceftazidime_avibactam': predict_cazavi_phenotype(blast_results, unified_results),
        'disclaimer': DISCLAIMER,
    }
