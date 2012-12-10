#!/bin/bash

while read -r line         
do     
    # the multiple -e are to work-around for MacOS BSD sed.
    # don't touch it unless you have a Mac to test it.
    pkg=$(echo ${line} | sed -n -e '/\\package_listing{[^}]*}/ { s/\\package_listing{\([^}]*\)}/\1/;' -e 'p' -e '}')

    if [ -n "${pkg}" ]; then
        top_level=${pkg%/*}
        lower_level=${pkg##*/}
        if [ -n "${top_level}" ]; then
            filename="../${top_level}/doc/${lower_level}/PackageDescription.txt"
        else
            filename="../${pkg}/doc/${pkg}/PackageDescription.txt"
        fi
        sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' < "$filename"
    else
        echo -E "${line}"
    fi
done <"$1"
