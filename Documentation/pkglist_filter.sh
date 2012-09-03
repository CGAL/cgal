#!/bin/sh

while read -r line         
do     
    pkg=$(echo ${line} | sed -n '/\\package_listing{[^}]*}/ { s/\\package_listing{\([^}]*\)}/\1/; p }')
    if [ -n "${pkg}" ]; then
        top_level=${pkg%/*}
        lower_level=${pkg##*/}
        if [ -n "${top_level}" ]; then
            filename="../${top_level}/doc/${lower_level}/PackageDescription.txt"
        else
            filename="../${pkg}/doc/${pkg}/PackageDescription.txt"
        fi
        sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' < "$filename" | \
            awk '/\\PkgDescriptionBegin{[^}]*}/ { match($0, "(\\\\PkgDescriptionBegin{)([^}]*)}", a); esc=a[2]; gsub(" ", "-", esc); printf("%s%s,%s}", a[1], a[2], esc); next} {print}'
    else
        echo -E "${line}"
    fi
done <"$1"
