all: normal cxx0x

CGAL_INC = -I/home/glisse/repos/cgal/trunk/Algebraic_foundations/include -I/home/glisse/repos/cgal/trunk/STL_Extension/include -I/data/repos/cgal/trunk/Algebraic_foundations/include -I/data/repos/cgal/trunk/STL_Extension/include -DCGAL_DISABLE_ROUNDING_MATH_CHECK

EIGEN_INC = `pkg-config --cflags eigen3|sed -e 's/-I/-isystem/g'` -DCGAL_USE_EIGEN

normal:
	g++ test.cpp -O2 -g -lCGAL -lboost_thread -frounding-math -Wall -I. -lmpfr -lgmp -I../../include ${EIGEN_INC} ${CGAL_INC} -DCGAL_CFG_NO_CPP0X_COPY_N -DCGAL_CFG_NO_CPP0X_NEXT_PREV
	./a.out

cxx0x:
	g++ -std=c++0x -O2 -g test.cpp -lCGAL -lboost_thread -frounding-math -Wall -I. -lmpfr -lgmp -DCGAL_CXX0X -I../../include ${EIGEN_INC} ${CGAL_INC}
	./a.out

#-DBOOST_RESULT_OF_USE_DECLTYPE 
