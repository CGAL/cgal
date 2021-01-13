#!/bin/bash
sudo add-apt-repository ppa:mikhailnov/pulseeffects -y
sudo apt-get update
sudo apt-get install -y  libmpfr-dev \
 libeigen3-dev qtbase5-dev libqt5sql5-sqlite libqt5opengl5-dev qtscript5-dev \
 libqt5svg5-dev qttools5-dev qttools5-dev-tools libboost1.72-dev zsh
#update cmake to 3.18.4
sudo apt purge --auto-remove cmake
cd /tmp
wget https://cmake.org/files/v3.18/cmake-3.18.4-Linux-x86_64.sh
sudo sh cmake-3.18.4-Linux-x86_64.sh --skip-license --prefix=/usr/local
rm cmake-3.18.4-Linux-x86_64.sh
