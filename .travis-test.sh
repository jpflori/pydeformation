#! /bin/sh
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e
SAGE_IMAGE=`python2 -c "import sage_version; print sage_version.get_all_version_names('${SAGE_SERVER}index.html',${SAGE_AGE})"`
$HOME/SageMath/sage -pip install --user --upgrade -v -i https://pypi.python.org/pypi sagemath # Check that Sage is installed
$HOME/SageMath/sage -pip install --upgrade --no-index -v .
$HOME/SageMath/sage setup.py test
(cd docs && $HOME/SageMath/sage -sh -c "make html")
$HOME/SageMath/sage -pip uninstall .
