#!/bin/bash

# This application is a cross-platform version of cgal_test.
# This is a script for the CGAL test suite. Such a script must obey
# the following rules:
#
# - for every target two one line messages are written to the file 'error.txt'
#     the first one indicates if the compilation was successful
#     the second one indicates if the execution was successful
#   if one of the two was not successful, the line should start with 'ERROR:'
# - running the script should not require any user interaction
# - applications must be already compiled

ERRORFILE=error.txt

#---------------------------------------------------------------------#
#                   find_executable <target>
#                   (different on Windows and Unix)
#---------------------------------------------------------------------#

find_executable()
{
    PARAM_APPLICATION=""
    [ -f ./debug/$1.exe ] && PARAM_APPLICATION="./debug/$1.exe"
    [ -f ./release/$1.exe ] && PARAM_APPLICATION="./release/$1.exe"
    [ -x ./$1 ] && PARAM_APPLICATION="./$1"
    echo "$PARAM_APPLICATION"
}

#---------------------------------------------------------------------#
#                    run <target>
#---------------------------------------------------------------------#

run()
{
    # Find exe
    COMMAND="`find_executable $1`"
    if [ -f "$COMMAND" ]; then
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
      ulimit -t 3600 2> /dev/null
      if eval $COMMAND 2>&1 ; then
          echo "   successful execution   of $1" >> $ERRORFILE
      else
          echo "   ERROR:     execution   of $1" >> $ERRORFILE
      fi
    else
      echo   "   ERROR:     not executed   $1" >> $ERRORFILE
    fi
}


#---------------------------------------------------------------------#
#                    main
#---------------------------------------------------------------------#

# start redirection to log file
(

#---------------------------------------------------------------------#
#                    remove the previous error file
#---------------------------------------------------------------------#

rm -f $ERRORFILE
touch $ERRORFILE

#---------------------------------------------------------------------#
#                    run the tests
#---------------------------------------------------------------------#

if [ $# -ne 0 ] ; then
  for file in $* ; do
    run $file
  done
else
  echo "Run all tests."
  run Authalic_parameterization 
  run Complete_parameterization_example 
  run Mesh_cutting_parameterization 
  run polyhedron_ex_parameterization
  run Simple_parameterization 
  run Square_border_parameterization 
fi

#---------------------------------------------------------------------#
#                   Recap results
#---------------------------------------------------------------------#

echo "------------------------------------------------------------------"
echo "- Results"
echo "------------------------------------------------------------------"
echo
cat $ERRORFILE
echo
rm -f $ERRORFILE

) 2>&1 | tee quick_test_suite.log

