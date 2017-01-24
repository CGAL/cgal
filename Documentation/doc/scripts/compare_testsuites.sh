#!/bin/bash
#Path to the CGAL_Documentation_build_directory/doc_output
PATH_TO_DOC="$1"

if ! [ -d "$PATH_TO_DOC" ] || [ $(basename $PATH_TO_DOC) != "doc_output" ]; then
  echo "wrong path"
  exit 1
fi

#path to the repository containing the output of this script for the reference documentation.
DOC_REF="$2"
#output in a new directory
mkdir -p doc_data
cd ./doc_data
if [ $# -gt 2 ]; then
  echo "too many arguments"
  exit 1
 fi

FAILURES=()
for dir in $PATH_TO_DOC/*
do
  OUTPUT=$(basename $dir)
  python ../documentation_parser.py $dir/xml > ./"$OUTPUT.txt"
  if [ $? -eq 0 ]; then
    echo "$dir OK"
  else
    echo "$dir FAILED"
    FAILURES+="$dir "
  fi
done
cd ..
echo "Output generated"
if ! [ -d "$DOC_REF" ]; then
  echo "No reference given. Script is finished."
  exit 2
fi

#diff the output and the reference output, ignoring the differences in whitespaces
diff -u -N -w ./doc_data $DOC_REF > ./diff.txt

#generate an html page showing the status of the diff
echo "<!DOCTYPE html>" > result.html
echo "<html>" >> result.html
echo "    <head>" >> result.html
echo "        <meta charset=\"utf-8\" />" >> result.html
echo "        <title>Documentation Status</title>" >> result.html
echo "    </head>" >> result.html
echo "    <body>" >> result.html
    
#if there is a diff, give a link to show it
if [[ -s ./diff.txt ]]
then
  echo "<p> The documentation has changed ! <a href=\"diff.txt\">See logs. </a> <br /> " >> result.html
#else just say that everything is fine
else
  echo "<p> The documentation has not changed. <br /><br /> " >> result.html
fi
if [ "${FAILURES[0]}" != "" ]; then
  echo " Some packages encountered problems while being parsed : <br /><br /> " >> result.html
  echo " ${FAILURE[*]} <br /> ">> result.html
fi
echo "    </p>" >> result.html
echo "    </body>" >> result.html
echo "</html>" >> result.html
