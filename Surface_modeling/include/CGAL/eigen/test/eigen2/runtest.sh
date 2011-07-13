#!/bin/bash

black='\E[30m'
red='\E[31m'
green='\E[32m'
yellow='\E[33m'
blue='\E[34m'
magenta='\E[35m'
cyan='\E[36m'
white='\E[37m'

if make test_$1 > /dev/null  2> .runtest.log ; then
  if ! ./test_$1 r20 > /dev/null 2> .runtest.log ; then
    echo -e  $red Test $1 failed: $black
    echo -e $blue
    cat .runtest.log
    echo -e $black
    exit 1
  else
  echo -e $green Test $1 passed$black
  fi
else
  echo -e  $red Build of target $1 failed: $black
  echo -e $blue
  cat .runtest.log
  echo -e $black
  exit 1
fi
