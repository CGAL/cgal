#! /bin/sh

rsync -av --chmod=u+rwX gauguin:.autocgalrc autocgalrc
rsync -av --chmod=u+rwX gauguin:cgal/launch_testsuite.bat .
rsync -Cvr --chmod=u+rwX --exclude .vs --exclude '*.cpp' --exclude '*vcxproj*' --exclude '*.sln' --exclude '*.sdf' --exclude include --exclude Makefile --exclude \*.cmake --exclude bin --exclude lib --exclude config --exclude CMakeFiles --exclude src gauguin:cgal/reference_platforms .
chmod -R u=rwX .
