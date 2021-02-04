#!/bin/bash

#use this script from inside the build directory of the Polyhedron demo

declare config="Release"

declare target_directory="CGAL_demo_with_dlls"
if [[ ! -d "$target_directory" ]]
then
  mkdir $target_directory
fi

copy_dll()
{
  local dll_full_path="`cygpath --unix --absolute $1`"
  echo "copy " $dll_full_path " to " $2
  cp $dll_full_path $2
}


files=($config/*.exe)
files+=($config/*.dll)
files+=(Plugins/*/$config/*.dll)

for file in "${files[@]}"; do

  # copy exe or dll
  copy_dll $file $target_directory

  # list and copy dependencies
  cygcheck $file | while read -r dll ; do
  
    copy_dll $dll $target_directory
  
  done; #check dependencies
    
done #loop over directories
