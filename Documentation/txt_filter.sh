#!/bin/sh

while read -r line         
do     
    pkg=$(echo ${line} | sed -n '/\\package_listing{[^}]*}/p' | sed 's/\\package_listing{\([^}]*\)}/\1/')
    if [ -n "${pkg}" ]; then
        cat "../${pkg}/doc/${pkg}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' | sed "s/\(\\PkgDescriptionBegin{[^}]*\)}/\1,${pkg}}/"
    else
        echo -E "${line}"
    fi
done <$1
