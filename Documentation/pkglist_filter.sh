#!/bin/sh

while read -r line         
do     
    pkg=$(echo ${line} | sed -n '/\\package_listing{[^}]*}/p' | sed 's/\\package_listing{\([^}]*\)}/\1/')
    if [ -n "${pkg}" ]; then
        top_level=${pkg%/*}
        lower_level=${pkg##*/}
        if [ -n "${top_level}" ]; then
            cat "../${top_level}/doc/${lower_level}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' | sed "s/\(\\PkgDescriptionBegin{[^}]*\)}/\1,${lower_level}}/"
        else
            cat "../${pkg}/doc/${pkg}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' | sed "s/\(\\PkgDescriptionBegin{[^}]*\)}/\1,${pkg}}/"
        fi
    else
        echo -E "${line}"
    fi
done <$1
