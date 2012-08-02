#!/bin/sh

while read -r line         
do     
    pkg=$(echo ${line} | sed -n '/\\package_listing{[^}]*}/p' | sed 's/\\package_listing{\([^}]*\)}/\1/')
    if [ -n "${pkg}" ]; then
        cat "../${pkg}/doc/${pkg}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p'
    else
        echo -E "${line}"
    fi
done <$1
