include ${CGAL_MAKEFILE}

OBJS=Gbrs_polynomial_1

all: ${OBJS}

Gbrs_polynomial_1: Gbrs_polynomial_1.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

clean:
	rm -f ${OBJS}

