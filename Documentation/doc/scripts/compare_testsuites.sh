#!/bin/bash
#Path to the CGAL_Documentation_build_directory/doc_output
PATH_TO_DOC="$1"
#path to the repository containing the output of this script for the reference documentation.
DOC_REF="$2"
#output in a new directory
mkdir -p doc_data
cd ./doc_data
for dir in $PATH_TO_DOC/*
do
  OUTPUT=$(basename $dir)
  python ../documentation_parser.py $dir/xml > ./"$OUTPUT.txt"
  echo "$dir done"
done
cd ..
#diff the output and the reference output, ignoring the differences in whitespaces
diff -u -w ./doc_data $DOC_REF > ./diff.txt

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
  echo " The documentation has changed ! <a href=\"diff.txt\">See logs. </a> " >> result.html
#else just say that everything is fine
else
  echo " The documentation has not changed. " >> result.html
fi
echo "    </body>" >> result.html
echo "</html>" >> result.html
