#! /bin/bash
#
# This script is to be run in COREBLAS's top-level directory.
# It will create a release file in the current directory.
#
# The script will search for the version numbers in the include directory.
# The script will check if the version numbers match in the CMake installation script.
# The script relies on using a modern TAR with advanced options that minimize the size of the resulting archive file.
#

MJR=`grep -r COREBLAS_VERSION_MAJOR include | awk '{print $NF;}'`
MNR=`grep -r COREBLAS_VERSION_MINOR include | awk '{print $NF;}'`
PTC=`grep -r COREBLAS_VERSION_PATCH include | awk '{print $NF;}'`

# this can come from include/coreblas_types.h
VERSION=${MJR}.${MNR}.$PTC

if test -z "`grep -i coreblas.version CMakeLists.txt | grep $VERSION`" ; then
  echo Version mismatch between headers $VERSION and CMakeLists.txt
  grep -i coreblas.version CMakeLists.txt
  exit 127
fi

DIR=coreblas-${VERSION}

if test ! -e $DIR ; then
    ln -s . $DIR
fi

echo Preparing $DIR ...

find -H ${DIR} -maxdepth 1 -type f -name '[A-Za-z0-9]*' | \
    xargs echo ${DIR}/*/*.hin ${DIR}/*/*.[hc] ${DIR}/config/*.py ${DIR}/tools/*.py ${DIR}/cmake/*.cmake ${DIR}/*/*.f90 ${DIR}/*/*.lua ${DIR}/*/doxygen* ${DIR}/lib/pkgconfig/coreblas.pc.in |  \
    xargs tar --exclude=.hgtags --exclude=coreblas_config.h --exclude=Makefile.\*.gen --owner=root --group=root --mtime=1970-01-01 -chof ${DIR}.tar
gzip --best --rsyncable --verbose ${DIR}.tar

rm -r $DIR
