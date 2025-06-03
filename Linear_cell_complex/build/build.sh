rm -rf debug
mkdir debug
cd debug && cmake -DCMAKE_OSX_ARCHITECTURES=arm64 -DCGAL_DIR=$CMAKE_INSTALLED_PREFIX/lib/CGAL -DCMAKE_BUILD_TYPE=Release -DCGAL_DATA_DIR=/Users/yamazakisouichirou/Documents/GSoC/cgal/Linear_cell_complex/examples/Linear_cell_complex/data -DCMAKE_PREFIX_PATH=/Users/yamazakisouichirou/Documents/GSoC/qt-everywhere-src-6.5.3/qtbase /Users/yamazakisouichirou/Documents/GSoC/cgal/Linear_cell_complex/examples/Linear_cell_complex && make
