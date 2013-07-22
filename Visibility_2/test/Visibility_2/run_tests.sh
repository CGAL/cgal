#!/bin/sh
# Author: Francisc Bungiu
# E-mail: fbungiu@gmail.com

for file in *
do
	if ! [[ -d $file ]]
	then
		if [[ -x $file ]]
		then
		    echo "Executing '$file'..."
		    ./$file
		fi
	fi
done