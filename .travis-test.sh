#! /bin/sh
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e
cd $HOME
sagemath/sage -pip install --user --upgrade -v -i https://pypi.python.org/pypi sagemath # Check that Sage is installed
sagemath/sage -pip install --upgrade --no-index -v .
sagemath/sage setup.py test
(cd docs && $HOME/sagemath/sage -sh -c "make html")
sagemath/sage -pip uninstall .
