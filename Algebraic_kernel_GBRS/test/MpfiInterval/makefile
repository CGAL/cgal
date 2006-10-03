include ${CGAL_MAKEFILE}

OBJS=Mpfi

all: ${OBJS}

Mpfi: Mpfi.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< /usr/lib/libmpfi.a -lmpfr -o $@

clean:
	rm -f ${OBJS}

