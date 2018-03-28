#This script must be called from the CGAL root.
set -e
set -x
while test $# -gt 0
do
    case "$1" in
        --help) echo "Usage: $0 <doxygen_exe_path> "
        echo " $0 must be called from the CGAL root directory. It will compile documentation for all packages using doxygen_exe_path, and deduce "
        echo "their dependencies. It will then compare them with the previous ones and output 1if the dependencies has changed, "
        echo "0 otherwise."
        exit 0
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
for pkg_path in $CGAL_ROOT/*
do
  pkg=$(basename $pkg_path)
  if [ -f $pkg_path/package_info/$pkg/dependencies ]; then
    mv $pkg_path/package_info/$pkg/dependencies $pkg_path/package_info/$pkg/dependencies.old
  fi
done

cmake -DCGAL_ENABLE_CHECK_HEADERS=TRUE -DDOXYGEN_EXECUTABLE="$DOX_PATH" -DCGAL_COPY_DEPENDENCIES=TRUE -DCMAKE_CXX_FLAGS="-std=c++11" ..
make -j$(nproc --all) packages_dependencies
echo " Checks finished"
for pkg_path in $CGAL_ROOT/*
do
  pkg=$(basename $pkg_path)
  if [ -f "$pkg_path/package_info/$pkg/dependencies" ]; then
    PKG_DIFF=$(diff -N -w  "$pkg_path/package_info/$pkg/dependencies.old" "$pkg_path/package_info/$pkg/dependencies" || true)
    if [ -n "$PKG_DIFF" ]; then
      HAS_DIFF=TRUE
      echo "Differences in $pkg: $PKG_DIFF"
    else
      echo "No differencies in $pkg dependencies."
    fi
    if [ -f $pkg_path/package_info/$pkg/dependencies.old ]; then
      rm $pkg_path/package_info/$pkg/dependencies.old
    fi
  fi
done
echo " Checks finished"
cd $CGAL_ROOT
rm -r dep_check_build
if [ -n "$HAS_DIFF" ]; then
  echo " You should run cmake with options CGAL_ENABLE_CHECK_HEADERS and CGAL_COPY_DEPENDENCIES ON, make the target packages_dependencies and commit the new dependencies files."
  exit 1
else
  echo "The dependencies are up to date."
  exit 0
fi
