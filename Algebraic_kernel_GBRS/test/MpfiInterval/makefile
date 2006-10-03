include ${CGAL_MAKEFILE}

OBJS=Mpfi

all: ${OBJS}

Mpfi: Mpfi.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

clean:
	rm -f ${OBJS}

