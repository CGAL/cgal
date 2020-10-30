#!/bin/bash
sudo add-apt-repository ppa:mikhailnov/pulseeffects -y
sudo apt-get update
sudo apt-get install -y  libmpfr-dev \
 libeigen3-dev qtbase5-dev libqt5sql5-sqlite libqt5opengl5-dev qtscript5-dev \
 libqt5svg5-dev qttools5-dev qttools5-dev-tools libboost1.72-dev
git clone --branch v3.18.4 --depth 1 git://github.com/Kitware/CMake.git
cd ./CMake
mkdir build
cd build
cmake ..
make -j2
sudo make install
cd ../..
rm -rf ./CMake
