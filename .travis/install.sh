#!/bin/bash

[ -n "$CGAL_DEBUG_TRAVIS" ] && set -x
DONE=0
sudo add-apt-repository ppa:mikhailnov/pulseeffects -y
sudo apt-get update

while [ $DONE = 0 ]
do
  DONE=1 && sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes install clang-10 zsh \
flex bison cmake graphviz libgmp-dev libmpfr-dev libmpfi-dev zlib1g-dev libeigen3-dev  \
qtbase5-dev libqt5sql5-sqlite libqt5opengl5-dev qtscript5-dev libqt5svg5-dev qttools5-dev qttools5-dev-tools qml-module-qtgraphicaleffects libopencv-dev mesa-common-dev libmetis-dev libglu1-mesa-dev \
libboost1.72-dev || DONE=0 && sudo apt-get update
done
exit 0

