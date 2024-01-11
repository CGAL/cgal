#!/bin/bash

#use this script from inside the build directory of the Polyhedron demo
#Needs the Qt6_DIR env variable set to <Qt6_ROOT>/lib/cmake/Qt6

#No config : in autotest_cgal we use NMake as generator
#If using MSVC Generator, declare config="Release"


declare config="$2"
declare target_directory="$1"

if [[ ! -d "$target_directory" ]]
then
  mkdir $target_directory
fi

copy_dll()
{
  local dll_full_path=$(cygpath --unix --absolute "$1")
  echo "copy $dll_full_path to $2"
  #remove all dlls from system and Visual
  if ! [[ "$dll_full_path" =~ "api-ms-win" ]] && ! [[ "$dll_full_path" =~ "system32" ]]; then
    cp "$dll_full_path" "$2"
  fi
}


files=($PWD/$config/*.exe)
files+=($PWD/$config/*.dll)
files+=($PWD/Plugins/*/$config/*.dll)

for file in "${files[@]}"; do

  # copy exe or dll
  copy_dll "$file" "$target_directory"

  # list and copy dependencies
  cygcheck "$file" | while read -r dll ; do

    copy_dll "$dll" "$target_directory"

  done; #check dependencies
done #loop over directories
mkdir -p "$target_directory/platforms"
cp "$Qt6_INSTALLATION_DIR/plugins/platforms/qwindows.dll" "$target_directory/platforms"
