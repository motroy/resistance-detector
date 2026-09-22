"""Run the pipeline over every synthetic scenario and compare with ground truth.

``create_test_genomes.py`` writes each genome together with the phenotype it is
built to produce.  This test regenerates those genomes and checks that the
pipeline reproduces the declared expectations exactly - the closest thing to an
end-to-end accuracy check that does not need real sequencing data.
"""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from fos_cazavi.acquired import BlastDetector
from fos_cazavi.phenotype import predict_phenotypes

from .conftest import requires_blast

pytestmark = requires_blast

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope='module')
def scenario_dir(tmp_path_factory):
    directory = tmp_path_factory.mktemp('scenarios')
    subprocess.run([sys.executable, str(ROOT / 'create_test_genomes.py'), str(directory)],
                   check=True, capture_output=True, cwd=ROOT)
    return directory


def load_expectations(scenario_dir):
    with open(scenario_dir / 'expected_results.tsv') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_scenarios_were_generated(scenario_dir):
    assert len(load_expectations(scenario_dir)) >= 10


# Parametrised from the committed expectations so the count never drifts out of
# sync with create_test_genomes.py.
_COMMITTED = ROOT / 'test_data' / 'expected_results.tsv'
_SCENARIO_NAMES = ([row['Genome'] for row in
                    csv.DictReader(open(_COMMITTED), delimiter='\t')]
                   if _COMMITTED.exists() else [])


@pytest.mark.parametrize('genome', _SCENARIO_NAMES)
def test_scenario_matches_expected_phenotype(genome, scenario_dir, database, tmp_path):
    expectations = {row['Genome']: row for row in load_expectations(scenario_dir)}
    expectation = expectations[genome]

    detector = BlastDetector(
        str(scenario_dir / expectation['Genome']), database,
        str(tmp_path / expectation['Genome']),
        organism=expectation['Organism'])
    results = detector.run()
    # Mirror the CLI exactly: the organism reaches the phenotype logic too.
    prediction = predict_phenotypes(results, organism=expectation['Organism'])

    assert prediction['fosfomycin']['phenotype'] == expectation['Expected_Fosfomycin'], (
        f"{expectation['Genome']}: {expectation['Rationale']}\n"
        f"evidence: {prediction['fosfomycin']['evidence']}")
    assert (prediction['ceftazidime_avibactam']['phenotype']
            == expectation['Expected_Ceftazidime_Avibactam']), (
        f"{expectation['Genome']}: {expectation['Rationale']}\n"
        f"evidence: {prediction['ceftazidime_avibactam']['evidence']}")
