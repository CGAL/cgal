#!/bin/bash
mkdir -p openmesh
cd openmesh
wget -O openmesh.tar.gz https://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.tar.gz
tar xf openmesh.tar.gz --strip-components=1
sed -i '94i #include <sys/time.h>' src/OpenMesh/Tools/Utils/conio.cc

mkdir build
cd build
cmake -DBUILD_APPS=FALSE ..
make -j2
sudo make -j2 install &>/dev/null

#clean up
cd ../..
rm -rf ./openmesh
