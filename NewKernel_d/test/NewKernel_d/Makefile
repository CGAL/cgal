all: normal cxx0x

CGAL_INC = -I/home/glisse/repos/cgal/trunk/Algebraic_foundations/include -I/home/glisse/repos/cgal/trunk/STL_Extension/include -I/data/repos/cgal/trunk/Algebraic_foundations/include -I/data/repos/cgal/trunk/STL_Extension/include

normal:
	g++ test.cpp -O2 -g -lCGAL -frounding-math -Wall -I. -lgmp -I../../include -I/usr/include/eigen2 -DCGAL_USE_EIGEN ${CGAL_INC}
	./a.out

cxx0x:
	/usr/lib/gcc-snapshot/bin/g++ -std=c++0x -O2 -g test.cpp -lCGAL -frounding-math -Wall -I. -lgmp -DCGAL_CXX0X -I../../include -I/usr/include/eigen2 -DCGAL_USE_EIGEN ${CGAL_INC}
	./a.out

#-DBOOST_RESULT_OF_USE_DECLTYPE 
