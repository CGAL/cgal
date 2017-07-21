#!/bin/bash
if [ "$1" == '--help' ]; then
  echo "Usage: $0 <doc 1> [doc 2]"
  exho "doc 1 and doc 2 are paths to doxygen outputs (doc_output)."
  echo "Parse the xml output of doc 1 and creates a directory with organized text files."
  echo "Then, if doc_2 is specified, do the same for its xml output and make the diff between them."
  exit 0
fi
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
if [ $# -gt 5 ]; then
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
  exit 0
fi

#diff the output and the reference output, ignoring the differences in whitespaces
diff -u -N -w ./doc_data $DOC_REF > ./diff.txt || true
