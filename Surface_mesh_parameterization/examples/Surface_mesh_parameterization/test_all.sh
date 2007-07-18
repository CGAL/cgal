#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder

./test_model.sh holes
./test_model.sh mannequin-devil
./test_model.sh mask_cone
./test_model.sh nefertiti
./test_model.sh rotor
./test_model.sh sphere966
./test_model.sh three_peaks

echo "***************************************************************"

for TST in Authalic_parameterization \
           Simple_parameterization \
           Taucs_parameterization \
           Mesh_cutting_parameterization \
           Square_border_parameterization \
           Complete_parameterization_example \
           polyhedron_ex_parameterization
do
    echo " "
    echo "*** $TST ***"

    # Find executable name (different on Windows and Unix)
    [ -f ./VC/debug/$TST.exe ] && PARAM_APPLICATION="./VC/debug/$TST.exe"
    [ -f ./VC/release/$TST.exe ] && PARAM_APPLICATION="./VC/release/$TST.exe"
    [ -f ./VC/x64/debug/$TST.exe ] && PARAM_APPLICATION="./VC/x64/debug/$TST.exe"
    [ -f ./VC/x64/release/$TST.exe ] && PARAM_APPLICATION="./VC/x64/release/$TST.exe"
    [ -x ./$TST ] && PARAM_APPLICATION="./$TST"
    if [ -z "$PARAM_APPLICATION" ]; then
        echo "Cannot find $TST executable"
        exit 1
    fi

    COMMAND="$PARAM_APPLICATION `cat $TST.cmd | tr '\n' ' ' | tr '\r' ' '`"
    eval $COMMAND 2>&1
done
