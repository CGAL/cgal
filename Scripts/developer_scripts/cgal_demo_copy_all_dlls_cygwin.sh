#!/bin/bash

#use this script from inside the build directory of the Polyhedron demo
#Needs the Qt5_DIR env variable set to <Qt5_ROOT>/lib/cmake/Qt5

#No config : in autotest_cgal we use NMake as generator
#If you are using Visual as Generator, declare config="Release"


declare config="$PWD"
declare target_directory="$1"

if [[ ! -d "$target_directory" ]]
then
  mkdir $target_directory
fi

copy_dll()
{
  local dll_full_path=$(cygpath --unix --absolute "$1")
  echo "copy $dll_full_path to $2"
  cp "$dll_full_path" "$2"
}


files=($config/*.exe)
files+=($config/*.dll)
files+=(Plugins/*/$config/*.dll)

for file in "${files[@]}"; do

  # copy exe or dll
  copy_dll "$file" "$target_directory"

  # list and copy dependencies
  cygcheck "$file" | while read -r dll ; do

    copy_dll "$dll" "$target_directory"

  done; #check dependencies
  mkdir -p "$target_directory/platforms"
  cp "$Qt5_DIR/../../../plugins/platforms/qwindows.dll" "$target_directory/platforms"
done #loop over directories
