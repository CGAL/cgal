include ${CGAL_MAKEFILE}

MYFLAGS=-g #-DCGAL_RS_DEBUG -DCGAL_RS_VERB=0

OBJS=Gbrs_polynomial_1 Gbrs_polynomial_2

all: ${OBJS}

Gbrs_polynomial_1: Gbrs_polynomial_1.C
	${CGAL_CXX} ${MYFLAGS} ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

Gbrs_polynomial_2: Gbrs_polynomial_2.C
	${CGAL_CXX} ${MYFLAGS} ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

clean:
	rm -f ${OBJS}

