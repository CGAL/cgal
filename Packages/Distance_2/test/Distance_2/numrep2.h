#if TESTR == 1
    typedef double testnt;
    typedef CGAL_Cartesian<testnt> TestR;
#endif

#if TESTR == 2
    typedef double testnt;
    typedef CGAL_Homogeneous<testnt> TestR;
#endif

#if TESTR == 3
    typedef CGAL_TestfieldC testnt;
    typedef CGAL_Cartesian<testnt> TestR;
#endif

#if TESTR == 4
    typedef CGAL_TestrepH testnt;
    typedef CGAL_Homogeneous<testnt> TestR;
#endif

#if TESTR == 5
#include <LEDA/REDEFINE_NAMES.h>
    typedef integer testnt;
#include <LEDA/UNDEFINE_NAMES.h>
    typedef CGAL_Homogeneous<testnt> TestR;
#endif
