all: normal cxx11

CXX = g++

CGAL_INC = -I. -I../../include -I../../../Algebraic_foundations/include -I../../../STL_Extension/include -I../../../Number_types/include -I../../../Kernel_23/include -I../../../Installation/include -DCGAL_DISABLE_ROUNDING_MATH_CHECK

EIGEN_INC = `pkg-config --cflags eigen3|sed -e 's/-I/-isystem/g'` -DCGAL_EIGEN3_ENABLED

normal:
	${CXX} Epick_d.cpp -O2 -lCGAL -lboost_thread -frounding-math -Wall -Wextra -lmpfr -lgmp ${CGAL_INC} ${EIGEN_INC}
	./a.out

cxx11:
	${CXX} -std=c++0x -O2 Epick_d.cpp -lCGAL -lboost_thread -frounding-math -Wall -Wextra -lmpfr -lgmp ${CGAL_INC} ${EIGEN_INC} -o b.out
	./b.out

#-DBOOST_RESULT_OF_USE_DECLTYPE 
