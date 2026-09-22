#!/usr/bin/env sh
set -eu

fos-cazavi batch -i .validation_genomes/PRJNA741867 -o bioproject_tests/PRJNA741867_test_results --organism Klebsiella_pneumoniae -j 6
fos-cazavi batch -i .validation_genomes/PRJNA595047 -o bioproject_tests/PRJNA595047_test --organism Klebsiella_pneumoniae -j 6
fos-cazavi batch -i .validation_genomes/PRJNA781811 -o bioproject_tests/PRJNA781811_test --organism Klebsiella_pneumoniae -j 6
fos-cazavi batch -i .validation_genomes/PRJNA1086695 -o bioproject_tests/PRJNA1086695_test --organism Klebsiella_pneumoniae -j 6
fos-cazavi batch -i .validation_genomes/Paeruginosa_ML_subset -o bioproject_tests/Paeruginosa_ML_subset --organism Pseudomonas_aeruginosa -j 6
fos-cazavi batch -i .validation_genomes/CREC_fosA3_China -o bioproject_tests/CREC_fosA3_China --organism Escherichia -j 6
fos-cazavi combine -i bioproject_tests/*_test bioproject_tests/*_test_results bioproject_tests/Paeruginosa_ML_subset bioproject_tests/CREC_fosA3_China -o bioproject_tests/all_validation
