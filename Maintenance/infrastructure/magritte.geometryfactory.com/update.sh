#! /bin/sh

rsync -av --chmod=u+rwX cgaltester@magritte:.autocgalrc autocgalrc
rsync -Cvr --chmod=u+rwX --exclude include --exclude Makefile --exclude \*.cmake --exclude bin --exclude lib --exclude config --exclude CMakeFiles --exclude src cgaltester@magritte:cgal_test/reference_platforms .
