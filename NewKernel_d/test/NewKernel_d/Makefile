all: normal cxx0x

CGAL_INC = -I../../../Algebraic_foundations/include -I../../../STL_Extension/include -I../../../Number_types/include -I../../../Kernel_23/include -DCGAL_DISABLE_ROUNDING_MATH_CHECK

EIGEN_INC = `pkg-config --cflags eigen3|sed -e 's/-I/-isystem/g'` -DCGAL_EIGEN3_ENABLED

normal:
	g++ test.cpp -O2 -lCGAL -lboost_thread -frounding-math -Wall -Wextra -I. -lmpfr -lgmp -I../../include ${EIGEN_INC} ${CGAL_INC} -DCGAL_CFG_NO_CPP0X_COPY_N -DCGAL_CFG_NO_CPP0X_NEXT_PREV -DCGAL_CFG_NO_CPP0X_EXPLICIT_CONVERSION_OPERATORS
	./a.out

cxx0x:
	g++ -std=c++0x -O2 test.cpp -lCGAL -lboost_thread -frounding-math -Wall -Wextra -I. -lmpfr -lgmp -DCGAL_CXX0X -I../../include ${EIGEN_INC} ${CGAL_INC} -o b.out
	./b.out

#-DBOOST_RESULT_OF_USE_DECLTYPE 
