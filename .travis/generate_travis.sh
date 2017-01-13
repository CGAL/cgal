#!/bin/bash
set -e
cd ../

if [ -f "$PWD/.travis/packages.txt" ] 
then
	rm "$PWD/.travis/packages.txt"
fi

#find all the packages
PACKAGES=()
for f in *
do
  if [ -d  "$f/examples/$f" ] || [ -d  "$f/test/$f" ] || [ -d  "$f/demo/$f" ]
	then
		PACKAGES+="$f "
  	echo "$f " >> ./.travis/packages.txt
	fi
done

if [ -f ".travis.yml" ] 
then
  #copy the current .travis.yml for later check
	cp ./.travis.yml ./.travis.old
  rm .travis.yml
fi
#writes the first part of the file	
old_IFS=$IFS  
IFS=$'\n'     
for LINE in $(cat "$PWD/.travis/template.txt")
do
	if [ "$LINE" != " matrix: " ]
	then
		echo "$LINE" >> .travis.yml
  else
  	break
	fi
done
echo " matrix: " >> .travis.yml
IFS=$' '
#writes the matrix
for package in ${PACKAGES[@]}
do
 	echo "  - PACKAGE='$package' " >> .travis.yml
done
 	echo "  - PACKAGE='CHECK' " >> .travis.yml
IFS=$'\n'     
#writes the end of the file
COPY=0
for LINE in $(cat "$PWD/.travis/template.txt")
do
	if [ "$LINE" == "install: " ]
	then
		COPY=1
	fi
	if [ $COPY == 1 ]
	then
		echo "$LINE" >> .travis.yml
	fi
done
IFS=$' '
#check if there are differences between the files
read -a DIFF <<<$(diff -q ./.travis.yml ./.travis.old)
if [ "${DIFF[0]}" != "" ]
then
  echo ".travis.yml has changed"
fi
#erase old travis
rm ./.travis.old
IFS=$old_IFS
