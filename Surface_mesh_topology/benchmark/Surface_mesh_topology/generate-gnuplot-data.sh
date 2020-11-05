#!/bin/bash

# Requirement:
# The three files res-quadrangulation-computation.txt, res-path-homotopy.txt, res-polygonal-schema.txt must exist (result of the script run-benchmarks.sh)

ficIn="res-path-homotopy.txt"

function extract_info_for_reduced_map()
{
    if [ $# -ne 2 ]
    then
        echo "ERROR in extract_info_for_reduced_map: need two arguments."
        exit 1
    fi
    
    ficIn="${1}"
    ficOut="${2}"

    # 0) Copy old file
    if [ -f "${ficOut}" ]
    then
        cp -f "${ficOut}" "${ficOut}.COPY" 
    fi
       
    # 1) Extract 2-map size and simplification times
    echo "# ==========#INITIAL-MAP==========     ==========#REDUCED-MAP========" > "${ficOut}" 
    echo "#darts   #vertices #edges   #faces     #darts #vertices #edges #faces    GlobalTime" >> "${ficOut}" 
    
    # 2.1) Extract the different info.
    grep "Initial map:" "${ficIn}" | cut -d '=' -f 2 | cut -d ',' -f 1 > restmp1.txt
    grep "Initial map:" "${ficIn}" | cut -d '=' -f 3 | cut -d ',' -f 1 > restmp2.txt
    grep "Initial map:" "${ficIn}" | cut -d '=' -f 4 | cut -d ',' -f 1 > restmp3.txt
    grep "Initial map:" "${ficIn}" | cut -d '=' -f 5 | cut -d ',' -f 1 > restmp4.txt
    
    grep "Reduced map:" "${ficIn}" | cut -d '=' -f 2 | cut -d ',' -f 1 > restmp5.txt
    grep "Reduced map:" "${ficIn}" | cut -d '=' -f 3 | cut -d ',' -f 1 > restmp6.txt
    grep "Reduced map:" "${ficIn}" | cut -d '=' -f 4 | cut -d ',' -f 1 > restmp7.txt
    grep "Reduced map:" "${ficIn}" | cut -d '=' -f 5 | cut -d ',' -f 1 > restmp8.txt
    
    grep "\[TIME\] Total time for computation of reduced map: " "${ficIn}" | cut -d ' ' -f 9  > restmp9.txt
    
    # 2.2) Regroup the different info in different columns of a same file.
    paste -d '\t' restmp1.txt restmp2.txt restmp3.txt restmp4.txt restmp5.txt restmp6.txt restmp7.txt restmp8.txt restmp9.txt >> "${ficOut}"

    # 3) Erase temp files
    rm -f restmp1.txt restmp2.txt restmp3.txt restmp4.txt restmp5.txt restmp6.txt restmp7.txt restmp8.txt restmp9.txt
}

function extract_time_path_homotopy()
{
    if [ $# -ne 2 ]
    then
        echo "ERROR in extract_time_path_homotopy: need two arguments."
        exit 1
    fi

    ficIn="${1}"
    ficOut="${2}"
    
    # 0) Copy old file
    if [ -f "${ficOut}" ]
    then
        cp -f "${ficOut}" "${ficOut}.COPY" 
    fi
    
    # 1) Extract computation times of homotopy tests.
    
    # 1.1) Extract the different info.
    echo "Path1" > restmp1.txt
    grep "Random seed: " "${ficIn}" | cut -d ' ' -f 6 >> restmp1.txt
    echo "Path2" > restmp2.txt
    grep "Random seed: " "${ficIn}" | cut -d ' ' -f 12 >> restmp2.txt
    echo "TimeContractible" > restmp3.txt
    grep "\[TIME\] is_contractible: " "${ficIn}" | cut -d ' ' -f 3 >> restmp3.txt
    echo "TimeHomotopy" > restmp4.txt
    grep "\[TIME\] are_freely_homotopic: " "${ficIn}" | cut -d ' ' -f 3 >> restmp4.txt
    
    # 1.2) Regroup the different info in different columns of a same file.
    echo "# Size_of_path1; Size_of_path 2; time_of_contractible_test; time_of_homotopy_test." > "${ficOut}"
    paste -d '\t' restmp1.txt restmp2.txt restmp3.txt restmp4.txt >> "${ficOut}"
    
    # 2) Erase temp files
    rm -f restmp1.txt restmp2.txt restmp3.txt restmp4.txt
}
 
# Bench 1
extract_info_for_reduced_map "res-quadrangulation-computation.txt" "computation-time-reduce-surface.dat"

# Bench 2
extract_time_path_homotopy "res-path-homotopy.txt" "computation-time-path-homotopy.dat"

# Bench 3
extract_time_path_homotopy "res-polygonal-schema.txt" "computation-time-polygonal-schema.dat"
