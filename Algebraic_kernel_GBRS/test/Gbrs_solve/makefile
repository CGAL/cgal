include ${CGAL_MAKEFILE}

OBJS=rs

all: ${OBJS}

rs: rs.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< ${CGAL_LDFLAGS} -o $@

clean:
	rm -f ${OBJS}

