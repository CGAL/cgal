#!/bin/bash

CHECK=
case $1 in
  --check) CHECK=y;;
esac

set -e
cd ../

if [ -f "$PWD/.travis/packages.txt" ] 
then
	rm "$PWD/.travis/packages.txt"
fi

#find all the packages
PACKAGES=()
INDEX=0
i=0
for f in *
do
  if [ -d  "$f/package_info/$f" ]
  then
    PACKAGES[$INDEX]+="$f "
	i=$[i+1]
	if [ $i = 3 ]
	then
	  i=0
 	  INDEX=$[INDEX+1]
	fi
  	echo "$f " >> ./.travis/packages.txt
  fi
done
if [ -f ".travis.yml" ] 
then
  #copy the current .travis.yml for later check
  mv ./.travis.yml ./.travis.old
fi
if [ -f ".appveyor.yml" ]
then
  mv ./.appveyor.yml ./.appveyor.old
fi
#writes the first part of the file	
old_IFS=$IFS
for TYPE in 'travis' 'appveyor'
do
IFS=$'\n'
  for LINE in $(cat "$PWD/.travis/$TYPE-template.txt")
  do
    if { [ $TYPE = 'travis' ] && [ "$LINE" != " include:" ]; } || \
{ [ $TYPE = 'appveyor' ] && [ "$LINE" != "  matrix:" ]; }
    then
      echo "$LINE" >> .$TYPE.yml
    else
      break
    fi
  done
  if [ $TYPE = 'appveyor' ]; then
    echo " matrix: " >> .$TYPE.yml
  else
    echo " include: " >> .$TYPE.yml
  fi
  #writes the matrix
  if [ $TYPE = 'travis' ]; then
    echo "  - compiler: gcc " >> .travis.yml
    echo "    env: PACKAGE='CHECK' " >> .travis.yml
    echo "  - compiler: clang-3.6" >> .travis.yml
    echo "    env: PACKAGE='CHECK' " >> .travis.yml
  fi
  for package in ${PACKAGES[@]}
  do
    if [ $TYPE = 'travis' ]; then
      echo "  - compiler: clang-3.6" >> .travis.yml
      echo "    env: PACKAGE='$package' " >> .travis.yml
    else
      echo "  - PACKAGE: '$package' " >> .$TYPE.yml
    fi
  done
  if [ $TYPE = 'travis' ]; then
    echo "  - compiler: clang-3.6" >> .travis.yml
    echo "    env: PACKAGE='Polyhedron_demo' " >> .travis.yml
  else
    echo "  - PACKAGE: 'Polyhedron_demo' " >> .$TYPE.yml
  fi

  #writes the end of the file
  COPY=0
  for LINE in $(cat "$PWD/.travis/$TYPE-template.txt")
  do
    if [ "$LINE" = "install:" ]
    then
      COPY=1
    fi
    if [ $COPY = 1 ]
    then
      echo "$LINE" >> .$TYPE.yml
    fi
  done
  IFS=$' '
  #check if there are differences between the files
  if ! cmp -s ./.$TYPE.yml ./.$TYPE.old;
  then
    echo ".$TYPE.yml has changed"
    if [ $TYPE = 'travis' ]; then
      if [ -n "$CHECK" ]; then
        echo "You should modify the file .travis/travis-template.txt"
        exit 1
      fi
    fi
  fi
  #erase old travis
  rm ./.$TYPE.old
IFS=$old_IFS
done

# Local Variables:
# tab-width: 2
# sh-basic-offset: 2
# End:
