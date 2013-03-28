#! /bin/sh

rsync -av --chmod=u+rwX picasso:.autocgalrc autocgalrc
rsync -av --chmod=u+rwX picasso:cgal/launch_testsuite.bat .
rsync -Cvr --chmod=u+rwX --exclude include --exclude Makefile --exclude \*.cmake --exclude bin --exclude lib --exclude config --exclude CMakeFiles --exclude src picasso:cgal/reference_platforms .
