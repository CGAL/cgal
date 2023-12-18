#!/bin/bash

set -ex

sudo apt-get update
sudo apt-get install -y \
 libmpfr-dev \
 libtbb-dev \
 libmetis-dev \
 libssh-dev \
 libeigen3-dev \
 qtbase5-dev libqt5sql5-sqlite libqt5opengl5-dev qtscript5-dev libqt5websockets5-dev \
 libqt5svg5-dev qttools5-dev qttools5-dev-tools \
 libboost-dev libboost-serialization-dev libboost-iostreams-dev libboost-filesystem-dev libboost-filesystem-dev \
 libvtk9-dev libgdcm-tools libvtkgdcm-dev libunwind-dev \
 libinsighttoolkit5-dev \
 libceres-dev \
 libglpk-dev \
 libopencv-dev \
 zsh \
 qt6-base-dev qt6-declarative-dev libqt6svg6-dev libqt6websockets6-dev

#update CMake
sudo apt purge --auto-remove cmake
cd /tmp
CMAKE_VER=$(curl --silent https://cmake.org/files/LatestRelease/cmake-latest-files-v1.json | jq -r .version.string)
wget https://cmake.org/files/LatestRelease/cmake-$CMAKE_VER-linux-x86_64.sh
sudo sh cmake-*.sh --skip-license --prefix=/usr/local
rm cmake-*.sh
