#if TESTR == 1
    typedef double testnt;
    typedef CGAL::Cartesian<testnt> TestR;
#endif

#if TESTR == 2
    typedef double testnt;
    typedef CGAL::Homogeneous<testnt> TestR;
#endif

#if TESTR == 3
    typedef CGAL::TestfieldC testnt;
    typedef CGAL::Cartesian<testnt> TestR;
#endif

#if TESTR == 4
    typedef CGAL::TestrepH testnt;
    typedef CGAL::Homogeneous<testnt> TestR;
#endif

#if TESTR == 5
#include <LEDA/REDEFINE_NAMES.h>
    typedef integer testnt;
#include <LEDA/UNDEFINE_NAMES.h>
    typedef CGAL::Homogeneous<testnt> TestR;
#endif
