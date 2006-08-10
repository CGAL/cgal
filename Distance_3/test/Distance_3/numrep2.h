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
#if CGAL_LEDA_VERSION < 500
#include <LEDA/REDEFINE_NAMES.h>
#else
#include <LEDA/internal/REDEFINE_NAMES.h>
#endif
    typedef integer testnt;
#if CGAL_LEDA_VERSION < 500
#include <LEDA/UNDEFINE_NAMES.h>
#else
#include <LEDA/internal/UNDEFINE_NAMES.h>
#endif
    typedef CGAL::Homogeneous<testnt> TestR;
#endif
