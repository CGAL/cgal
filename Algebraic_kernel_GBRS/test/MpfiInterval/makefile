include ${CGAL_MAKEFILE}

OBJS=test_int test_gmp test_mpfr test

all: ${OBJS}

test_int: test_int.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< /usr/lib/libmpfi.a -lmpfr -o $@

test_gmp: test_gmp.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< /usr/lib/libmpfi.a -lmpfr -o $@

test_mpfr: test_mpfr.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< /usr/lib/libmpfi.a -lmpfr -o $@

test: test.C
	${CGAL_CXX} -g ${CGAL_CXXFLAGS} ${CGAL_LIBPATHFLAGS} $< /usr/lib/libmpfi.a -lmpfr -o $@

clean:
	rm -f ${OBJS}

