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


@pytest.mark.parametrize('index', range(14))
def test_scenario_matches_expected_phenotype(index, scenario_dir, database, tmp_path):
    expectations = load_expectations(scenario_dir)
    if index >= len(expectations):
        pytest.skip('scenario not defined')
    expectation = expectations[index]

    detector = BlastDetector(
        str(scenario_dir / expectation['Genome']), database,
        str(tmp_path / expectation['Genome']),
        organism=expectation['Organism'])
    results = detector.run()
    prediction = predict_phenotypes(results)

    assert prediction['fosfomycin']['phenotype'] == expectation['Expected_Fosfomycin'], (
        f"{expectation['Genome']}: {expectation['Rationale']}\n"
        f"evidence: {prediction['fosfomycin']['evidence']}")
    assert (prediction['ceftazidime_avibactam']['phenotype']
            == expectation['Expected_Ceftazidime_Avibactam']), (
        f"{expectation['Genome']}: {expectation['Rationale']}\n"
        f"evidence: {prediction['ceftazidime_avibactam']['evidence']}")
