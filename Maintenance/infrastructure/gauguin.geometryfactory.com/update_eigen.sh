#!/bin/bash
cd /cygdrive/c/3rdPartyLibs/eigen-master
git pull origin master
cd ../boost_master
git submodule update --init --recursive
./b2 -j8 --toolset=msvc-14.1 address-model=64 architecture=x86 link=static --prefix="C:\3rdPartyLibs\boost_master\install_dir" install