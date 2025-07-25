#!/bin/sh
# creates an archive (tar.gz) of the source

NAME="StraightSkel"
VERSION=$(date +%Y%m%d)
TMP=${TMP:-/tmp}
CWD=$(pwd)

set -e   # exit on error

mkdir "${TMP}/${NAME}"
cp -r src test res doc *.in *.txt *.sh *.py "${TMP}/${NAME}"
rm -rf $(find "${TMP}/${NAME}" -type d -name '.svn')
cd "${TMP}"
tar czvf "${NAME}-${VERSION}.tar.gz" "${NAME}"
mv "${NAME}-${VERSION}.tar.gz" "${CWD}"
rm -rf "${NAME}"
