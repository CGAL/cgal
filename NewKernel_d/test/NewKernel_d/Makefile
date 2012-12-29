all: normal cxx0x

CGAL_INC = -I/home/glisse/repos/cgal/master/Algebraic_foundations/include -I/home/glisse/repos/cgal/master/STL_Extension/include -I/home/glisse/repos/cgal/master/Number_types/include -I/data/repos/cgal/master/Algebraic_foundations/include -I/data/repos/cgal/master/STL_Extension/include -I/data/repos/cgal/master/Number_types/include -DCGAL_DISABLE_ROUNDING_MATH_CHECK

EIGEN_INC = `pkg-config --cflags eigen3|sed -e 's/-I/-isystem/g'` -DCGAL_EIGEN3_ENABLED

normal:
	g++ test.cpp -O2 -lCGAL -lboost_thread -frounding-math -Wall -Wextra -I. -lmpfr -lgmp -I../../include ${EIGEN_INC} ${CGAL_INC} -DCGAL_CFG_NO_CPP0X_COPY_N -DCGAL_CFG_NO_CPP0X_NEXT_PREV
	./a.out

cxx0x:
	g++ -std=c++0x -O2 test.cpp -lCGAL -lboost_thread -frounding-math -Wall -Wextra -I. -lmpfr -lgmp -DCGAL_CXX0X -I../../include ${EIGEN_INC} ${CGAL_INC}
	./a.out

#-DBOOST_RESULT_OF_USE_DECLTYPE 
