#!/bin/sh

while read -r line         
do     
    pkg=$(echo ${line} | sed -n '/\\package_listing{[^}]*}/p' | sed 's/\\package_listing{\([^}]*\)}/\1/')
    if [ -n "${pkg}" ]; then
        top_level=${pkg%/*}
        lower_level=${pkg##*/}
        if [ -n "${top_level}" ]; then
            cat "../${top_level}/doc/${lower_level}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' | \
                awk '/\\PkgDescriptionBegin{[^}]*}/ { match($0, "(\\\\PkgDescriptionBegin{)([^}]*)}", a); esc=a[2]; gsub(" ", "-", esc); printf("%s%s,%s}", a[1], a[2], esc); next} {print}'
        else
            cat "../${pkg}/doc/${pkg}/PackageDescription.txt" | sed -n '/PkgDescriptionBegin/,/PkgDescriptionEnd/p' | \
                awk '/\\PkgDescriptionBegin{[^}]*}/ { match($0, "(\\\\PkgDescriptionBegin{)([^}]*)}", a); esc=a[2]; gsub(" ", "-", esc); printf("%s%s,%s}", a[1], a[2], esc); next} {print}'
        fi
    else
        echo -E "${line}"
    fi
done <$1
