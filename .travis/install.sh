#!/bin/bash

[ -n "$CGAL_DEBUG_TRAVIS" ] && set -x
DONE=0
while [ $DONE = 0 ]
do
  DONE=1 && sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes install clang zsh \
flex bison cmake graphviz libgmp-dev libmpfr-dev libmpfi-dev zlib1g-dev libeigen3-dev  libboost-dev \
libboost-system-dev libboost-program-options-dev libboost-thread-dev libboost-iostreams-dev \
qtbase5-dev qtscript5-dev libqt5svg5-dev qttools5-dev qttools5-dev-tools qml-module-qtgraphicaleffects libopencv-dev mesa-common-dev libmetis-dev libglu1-mesa-dev \
|| DONE=0 && sudo apt-get update
done
exit 0

