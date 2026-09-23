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
    AVIBACTAM_COMBINATION_SCOPE, CAZAVI_CONTRIBUTORY_GENES, CAZAVI_SCOPE,
    FOSA_FAMILIES, FOS_SCOPE, FOS_TRANSPORT_GENES, FOS_UNDEFINED_BREAKPOINT_ORGANISMS,
    INTRINSIC_FOSA_ORGANISMS, INTRINSIC_GENES,
    PERMEABILITY_AMPLIFIED_FAMILIES, PORIN_GENES, UNASSESSED_CAZAVI_MECHANISM,
    gene_family, is_mbl, mutation_scope,
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


def predict_fos_phenotype(blast_results, unified_results=None, organism=None):
    """Predict fosfomycin susceptibility."""
    resistant, uncertain = [], []

    # Klebsiella and P. aeruginosa always carry a chromosomal fosA.  If one was
    # recognised, any *additional* fosA-family gene is genuinely acquired.  If
    # none was, a lone fosA hit cannot be told apart from a divergent copy of
    # that chromosomal gene by sequence identity alone.
    intrinsic_found = any(result['gene'] in INTRINSIC_GENES
                          for result in blast_results or [])
    intrinsic_expected = organism in INTRINSIC_FOSA_ORGANISMS

    for result in blast_results or []:
        gene = result['gene']
        family = gene_family(gene)

        if gene in INTRINSIC_GENES:
            # Intrinsic chromosomal fosA (e.g. fosAKP in K. pneumoniae) is
            # present in susceptible isolates and is not scored here.
            continue

        if family in FOSA_FAMILIES:
            if intrinsic_expected and not intrinsic_found:
                uncertain.append(
                    f"fosA-family enzyme {result['allele']} detected "
                    f"({result['identity']}% identity), but no intrinsic "
                    f"chromosomal fosA was recognised in this {organism.replace('_', ' ')} "
                    f"genome. Every isolate of this species carries one, so this may "
                    f"be a divergent chromosomal enzyme rather than an acquired gene")
            else:
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
            else:
                # Only mutations curated *for fosfomycin* count here.  Several
                # of these genes also carry mutations curated for other drugs -
                # cyaA_S352T is fosmidomycin, galU_R101C is cephalosporin - and
                # scoring those would be a plain false positive.
                fosfomycin_mutations = [
                    row['Label'] for row in result.get('curated_mutations', [])
                    if mutation_scope(row) == FOS_SCOPE
                ]
                if fosfomycin_mutations:
                    resistant.append(
                        f"Curated fosfomycin-resistance mutation(s) in {gene}: "
                        f"{', '.join(fosfomycin_mutations)}")

    if organism in FOS_UNDEFINED_BREAKPOINT_ORGANISMS and not resistant:
        # No validated clinical breakpoint exists for this organism/drug pair
        # (EUCAST publishes only an ECOFF for Pseudomonas, explicitly not a
        # clinical breakpoint; CLSI does not cover it at all). Without a
        # concrete mechanism in hand, neither Susceptible nor Resistant is
        # supportable from genotype alone - asserting either would overstate
        # the evidence, so this is Indeterminate rather than a guess.
        species = organism.replace('_', ' ')
        uncertain.append(
            f"No validated clinical breakpoint for fosfomycin exists for "
            f"{species} (EUCAST publishes only an epidemiological cut-off, "
            f"not a clinical breakpoint, for Pseudomonas spp.); genotype alone "
            f"cannot support a Susceptible or Resistant call for this species")

    return _resolve(
        resistant, uncertain,
        'No acquired fosfomycin-modifying enzyme and no loss-of-function change '
        'in the fosfomycin uptake/regulatory genes were detected')


def _contributory_cazavi_evidence(blast_results, amplifiable_enzymes):
    """Chromosomal changes that raise ceftazidime-avibactam MICs without being
    sufficient on their own.

    OmpK36 loss, PBP3 (FtsI) changes and EnvZ changes act by reducing drug entry
    or altering the target, which amplifies a beta-lactamase rather than
    defeating avibactam. They are reported as uncertain evidence, so on their own
    they give Indeterminate, and alongside a KPC escape variant they appear as
    supporting context.
    """
    evidence = []

    for result in blast_results or []:
        gene = result['gene']
        if gene not in CAZAVI_CONTRIBUTORY_GENES:
            continue

        for row in result.get('curated_mutations', []):
            scope = mutation_scope(row)
            if scope not in (CAZAVI_SCOPE, AVIBACTAM_COMBINATION_SCOPE):
                continue
            subclass = (row.get('Subclass') or '').upper()
            evidence.append(
                f"{row['Label']} in {gene}: curated by AMRFinderPlus for "
                f"{subclass.replace('/', ', ')}. This contributes to raised "
                f"ceftazidime-avibactam MICs but is not on its own established "
                f"as conferring resistance")

        # A knocked-out porin is stronger evidence than any single substitution
        # in it, but only where it is documented to matter: alongside a
        # beta-lactamase whose activity reduced drug entry amplifies.
        if (gene in PORIN_GENES and amplifiable_enzymes
                and result.get('loss_of_function') and result.get('complete')):
            evidence.append(
                f"Loss of function in porin {gene} ({result['lof_description']}) "
                f"alongside {', '.join(sorted(amplifiable_enzymes))}: reduced drug "
                f"entry combined with beta-lactamase activity is a documented route "
                f"to raised ceftazidime-avibactam MICs, though it is not on its own "
                f"established as conferring resistance")

    return evidence


def predict_cazavi_phenotype(blast_results, unified_results=None, organism=None):
    """Predict ceftazidime-avibactam susceptibility."""
    resistant, uncertain = [], []
    kpc_seen = False
    # OmpK36 loss is scored only alongside a KPC, the context the literature
    # documents; it is computed up front because the porin hit can come first.
    amplifiable_enzymes = {
        result['gene'] for result in blast_results or []
        if gene_family(result['gene']) in PERMEABILITY_AMPLIFIED_FAMILIES
    }

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

    uncertain.extend(_contributory_cazavi_evidence(blast_results, amplifiable_enzymes))

    unassessed = UNASSESSED_CAZAVI_MECHANISM.get(organism)
    if unassessed and not resistant:
        # Without the dominant mechanism in view, "susceptible" would overstate
        # what was actually ruled out.
        uncertain.append(
            f"No mechanism this tool assesses was found, but {unassessed}")

    if kpc_seen and not resistant and not uncertain:
        nothing_found = ('blaKPC detected but carrying no Omega-loop, 237-243 or '
                         'insertion-loop change associated with avibactam escape, '
                         'and no contributory porin/PBP3/EnvZ change; avibactam is '
                         'expected to inhibit it')
    else:
        nothing_found = ('No ceftazidime-avibactam resistance mechanism detected '
                         '(no metallo-beta-lactamase, no blaKPC escape variant, no '
                         'contributory porin/PBP3/EnvZ change)')

    return _resolve(resistant, uncertain, nothing_found)


def predict_phenotypes(blast_results, unified_results=None, organism=None):
    return {
        'fosfomycin': predict_fos_phenotype(blast_results, unified_results, organism),
        'ceftazidime_avibactam': predict_cazavi_phenotype(
            blast_results, unified_results, organism),
        'disclaimer': DISCLAIMER,
    }
