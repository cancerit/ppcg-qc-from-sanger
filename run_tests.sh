#!/usr/bin/env bash
set -e
pytest --cov-branch --cov-report term --cov-report html --cov=ppcg_qc_from_sanger --cov-fail-under=67
set +e

# these should not die:

echo -e "\n#################################"
echo      "# Running pycodestyle (style)   #"
echo      "#################################"
pycodestyle ppcg_qc_from_sanger

echo -e "\n#########################################"
echo      "# Running radon (cyclomatic complexity) #"
echo      "#########################################"
radon cc -nc ppcg_qc_from_sanger

echo -e "\n#########################################"
echo      "# Running radon (maintainability index) #"
echo      "#########################################"
radon mi -s ppcg_qc_from_sanger

echo -e "\n##############################"
echo      "# Running mdl (markdownlint) #"
echo      "##############################"
mdl -r ~MD013 .  # ignore line length rule.

exit 0 # don't die based on assements of code quality
