#if TESTR == 1

typedef double testnt;
typedef CGAL::Cartesian<testnt> TestR;

#endif

#if TESTR == 2

typedef long testnt;
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

#include <CGAL/MP_Float.h>
//typedef double testnt; // some tests trigger a degenerate case
typedef CGAL::MP_Float testnt;
typedef CGAL::Homogeneous<testnt> TestR;

#endif

#if TESTR == 6

typedef CGAL::Checked_long testnt;
typedef CGAL::Homogeneous<testnt> TestR;

#endif
