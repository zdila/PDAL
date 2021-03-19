#!/bin/bash

conda install -c conda-forge conda-build anaconda-client ceres-solver -y
conda install -c saedrna geographiclib
pwd
ls
git clone https://github.com/conda-forge/pdal-feedstock.git

cd pdal-feedstock
cat > recipe/recipe_clobber.yaml <<EOL
source:
  path: ../../
  url:
  sha256:

build:
  number: 2112
EOL

ls recipe
