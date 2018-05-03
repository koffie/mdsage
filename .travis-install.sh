#! /bin/bash
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e

if [ "${SAGE_AGE}" == "-1" ]; then
 sudo add-apt-repository ppa:aims/sagemath -y
 sudo apt-get update -qq
 sudo apt-get install sagemath-upstream-binary -y
 cd $HOME
 mkdir -p SageMath
else
  SAGE_IMAGE=`python2 -c "import setup; print setup.get_all_version_names('${SAGE_SERVER}index.html',${SAGE_AGE},broken = ${SAGE_BROKEN})"`
  cd $HOME
  echo "Obtaining Sage image:" ${SAGE_IMAGE}
  if [ ! -x SageMath/sage ] ; then
      rm -f SageMath.tar.bz2
      wget --progress=dot:giga ${SAGE_SERVER}${SAGE_IMAGE} -O SageMath.tar.bz2
      tar xf SageMath.tar.bz2
  fi
  MAKE="make -j4"
  export MAKE
  # Install packages
  #SageMath/sage -i lrslib
  # To initialize matplotlib font manager
  $HOME/SageMath/sage -python -c 'import matplotlib.pyplot'
fi

