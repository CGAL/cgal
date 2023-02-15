#!/bin/bash
if [ "$1" == '--help' -o ! -d "$1" -o ! -d "$2" ]; then
  echo "Usage: $0 <path to CGAL release> <path to output> <number of cores to dedicate>"
  echo "Builds and packages the Polyhedron demo form the CGAL dir."
  exit 0
fi
docker run --rm -v "$2":/results:Z -v "$1":/cgal:ro,z -e "NUMBER_OF_DEDICATED_CORES=$3" docker.io/cgal/bundle-3d-demo

