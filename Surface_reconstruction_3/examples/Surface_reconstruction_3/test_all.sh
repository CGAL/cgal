#!/bin/bash

############################################################
# This application is a cross-platform version of cgal_test.
# Applications must be already compiled.
############################################################

ERRORFILE=error.txt

# Find executable name (different on Windows and Unix)
# ----------------------------------------------------
find_executable()
{
    PARAM_APPLICATION=""
    [ -f ./VC/debug/$1.exe ] && PARAM_APPLICATION="./VC/debug/$1.exe"
    [ -f ./VC/release/$1.exe ] && PARAM_APPLICATION="./VC/release/$1.exe"
    [ -f ./VC/x64/debug/$1.exe ] && PARAM_APPLICATION="./VC/x64/debug/$1.exe"
    [ -f ./VC/x64/release/$1.exe ] && PARAM_APPLICATION="./VC/x64/release/$1.exe"
    [ -x ./$1 ] && PARAM_APPLICATION="./$1"
    echo "$PARAM_APPLICATION"
}

# run 1 test
# ----------
run()
{
    # Find exe
    COMMAND="`find_executable $1`"
    if [ -z "$COMMAND" ]; then
        echo "Cannot find $1 executable"
        exit 1
    fi

    # Add params
    if [ -f $1.cmd ] ; then
        COMMAND="$COMMAND `cat $1.cmd`"
    fi
    if [ -f $1.cin ] ; then
        COMMAND="cat $1.cin | $COMMAND"
    fi

    # Run
    echo "------------------------------------------------------------------"
    echo "- Executing $1"
    echo "------------------------------------------------------------------"
    echo
    if eval $COMMAND 2>&1 ; then
        echo "   successful execution   of $1" >> $ERRORFILE
    else
        echo "   ERROR:     execution   of $1" >> $ERRORFILE
    fi
    echo
}

# main
# ----

(

# remove the previous error file
rm -f $ERRORFILE
touch $ERRORFILE

# run the tests
if [ $# -ne 0 ] ; then
  for file in $* ; do
    run $file
  done
else
  run APSS_reconstruction
  run normal_estimation
  run poisson_reconstruction
fi

# Recap results
echo "------------------------------------------------------------------"
echo "- Results"
echo "------------------------------------------------------------------"
echo
cat $ERRORFILE
echo
rm -f $ERRORFILE

) 2>&1 | tee test_all.log

