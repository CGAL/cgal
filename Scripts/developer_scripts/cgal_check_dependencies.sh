#!/bin/bash
#This script must be called from the CGAL root.
set -e
if [ -n "$RUNNER_DEBUG" ]; then
  set -x
  CMAKE_DEBUG_OPT=" -v"
fi

while test $# -gt 0
do
    case "$1" in
        --help) echo "Usage: $0 <doxygen_exe_path> "
        echo " $0 must be called from the CGAL root directory. It will compile documentation for all packages using doxygen_exe_path, and deduce "
        echo "their dependencies. It will then compare them with the previous ones and output 1if the dependencies has changed, "
        echo "0 otherwise."
        exit 0
            ;;
        --check_headers) DO_CHECK_HEADERS="True"
            ;;
        --*) echo "bad option $1"
            ;;
        *) DOX_PATH="$1"
            ;;
    esac
    shift
done

CGAL_ROOT=$PWD
mkdir -p dep_check_build && cd dep_check_build
for pkg_path in "$CGAL_ROOT"/*
do
  pkg=$(basename "$pkg_path")
  if [ -f "$pkg_path/package_info/$pkg/dependencies" ]; then
    mv "$pkg_path/package_info/$pkg/dependencies" "$pkg_path/package_info/$pkg/dependencies.old"
  else
    if [ -d "$pkg_path/package_info/$pkg" ]; then
      touch "$pkg_path/package_info/$pkg/dependencies.old"
    fi
  fi
done

group()
{
    local title=$1
    shift
    if [ "$GITHUB_ACTIONS" = "true" ]; then
        echo "::group::$title"
    else
        printf '\n## %s\n\n' "$title"
    fi
    local status=0
    "$@" || status=$?
    if [ "$GITHUB_ACTIONS" = "true" ]; then
        echo "::endgroup::"
    fi
    return "$status"
}

group "Configure dependency check" cmake -DCGAL_ENABLE_CHECK_HEADERS=TRUE -DDOXYGEN_EXECUTABLE="$DOX_PATH" -DCGAL_COPY_DEPENDENCIES=TRUE -DCMAKE_CXX_FLAGS="-std=c++1y" ..
if [ -n "$DO_CHECK_HEADERS" ]; then
    group "Check headers" cmake --build . -j"$(nproc --all)" --target check_headers ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
    group "Check headers linked twice" cmake --build . -j"$(nproc --all)" --target check_headers_linked_twice ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
fi
group "Collect package dependencies" cmake --build . -j"$(nproc --all)" --target packages_dependencies ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
echo " Checks finished"
for pkg_path in "$CGAL_ROOT"/*
do
  pkg=$(basename "$pkg_path")
  if [ -f "$pkg_path/package_info/$pkg/dependencies" ]; then
    DIFF_STATUS=0
    PKG_DIFF=$(diff -u "$pkg_path/package_info/$pkg/dependencies.old" "$pkg_path/package_info/$pkg/dependencies") || DIFF_STATUS=$?
    if [ "$DIFF_STATUS" -eq 1 ]; then
      TOTAL_RES="Differences in $pkg:\n$PKG_DIFF\n$TOTAL_RES"
    elif [ "$DIFF_STATUS" -ne 0 ]; then
      echo "Failed to compare dependencies for $pkg" >&2
      exit "$DIFF_STATUS"
    fi
    if [ -f "$pkg_path/package_info/$pkg/dependencies.old" ]; then
      rm "$pkg_path/package_info/$pkg/dependencies.old"
    fi
  fi
done
echo " Checks finished"
cd "$CGAL_ROOT"
rm -r dep_check_build
if [ -n "$TOTAL_RES" ]; then
  # shellcheck disable=SC2059
  printf "$TOTAL_RES"
  echo " You can run cmake with options \`CGAL_ENABLE_CHECK_HEADERS\` and \`CGAL_COPY_DEPENDENCIES\` set to \`ON\`"
  echo " then build the target \`packages_dependencies\` and commit the new dependencies files,"
  echo " or simply manually edit the problematic files."
  exit 1
else
  echo "The dependencies are up to date."
  exit 0
fi
