#!/bin/bash

set -x
DONE=0
while [ $DONE = 0 ]
do
  DONE=1 && sudo -E apt-add-repository -y "ppa:ppsspp/cmake" ||  DONE=0 && sleep 5
done
DONE=0
while [ $DONE = 0 ]
do
  DONE=1 && sudo -E apt-add-repository -y "ppa:hedges/qt5.5" ||  DONE=0 && sleep 5
done

for pkg in clang-3.6 zsh cmake libboost1.55-dev libboost-system1.55-dev libboost-program-options1.55-dev libboost-thread1.55-dev libboost-iostreams1.55-dev libgmp-dev libmpfr-dev libmpfi-dev zlib1g-dev libeigen3-dev qt55base qt55script qt55svg qt55tools qt55graphicaleffects mesa-common-dev libglu1-mesa-dev
do
  DONE=0
  while [ $DONE = 0 ]
  do
    DONE=1 && sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes install $pkg || DONE=0 && sudo apt-get update
  done
done
