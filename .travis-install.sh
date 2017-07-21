#! /bin/sh
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e
SAGE_IMAGE=`python2 -c "import sage_version; print sage_version.get_all_version_names('${SAGE_SERVER}index.html',${SAGE_AGE})"`
cd $HOME
if [ ! -x sagemath/sage ] ; then 
    rm -f sagemath.tar.bz2
    wget ${SAGE_SERVER}${SAGE_IMAGE} -O sagemath.tar.bz2
    tar xf sagemath.tar.bz2
fi
MAKE="make -j4"
export MAKE
# Install packages
sagemath/sage -i deformation
# To initialize matplotlib font manager
sagemath/sage -python -c 'import matplotlib.pyplot'
