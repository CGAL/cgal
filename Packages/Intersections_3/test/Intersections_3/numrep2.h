#if TESTR == 1

typedef double testnt;
typedef CGAL_Cartesian<testnt> TestR;

#endif

#if TESTR == 2

typedef long testnt;
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

typedef double testnt;
typedef CGAL_Homogeneous<testnt> TestR;

#endif

#if TESTR == 6

typedef CGAL_Checked_long testnt;
typedef CGAL_Homogeneous<testnt> TestR;

#endif
