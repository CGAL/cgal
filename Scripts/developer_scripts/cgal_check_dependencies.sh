#This script must be called from the CGAL root.

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
for pkg in ../*
do
  if [ -f $pkg/dependencies ]; then
    echo "$pkg"
    mv $pkg/dependencies $pkg/dependencies.old
  fi
done

cmake -DCGAL_ENABLE_CHECK_HEADERS=TRUE -DDOXYGEN_EXECUTABLE="$DOX_PATH" -DCGAL_COPY_DEPENDENCIES=TRUE ..
make -j$(nproc --all) packages_dependencies
for pkg in ../*
do
  if [ -f $pkg/dependencies ]; then
    PKG_DIFF=$(diff -N -w  $pkg/dependencies $pkg/dependencies.old)
    if [ -n "$PKG_DIFF" ]; then
      HAS_DIFF=TRUE
      echo "Differences in $pkg: $PKG_DIFF"
      fi
    rm $pkg/dependencies.old
  fi
done
cd $CGAL_ROOT
rm -r dep_check_build
if [ -n "$HAS_DIFF" ]; then
  echo " You should run cmake with option CGAL_CHECK_HEADERS ON, make the target packages_dependencies and commit the new dependencies files."
  exit 1
else
  echo "The dependencies are up to date."
  exit 0
fi
