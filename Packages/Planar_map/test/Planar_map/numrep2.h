#if TESTR == 1

typedef double testnt;
typedef CGAL::Simple_Cartesian<testnt>          Simple_TestR;
typedef CGAL::Filtered_kernel<Cartesian_TestR>  TestR;
typedef TestR::RT inputt;

#endif

#if TESTR == 2

//typedef double testnt; // not used any more due to robustness while using homogeneous.
//typedef leda_debug_integer testnt;
typedef leda_integer testnt;
typedef CGAL::Homogeneous<testnt> TestR;
typedef TestR::RT inputt;

#endif


#if TESTR == 3

typedef CGAL::Gmpz testnt;
typedef CGAL::Cartesian<testnt> TestR;
typedef TestR::RT inputt;
#endif


#if TESTR == 4

typedef CGAL::Gmpz testnt;
typedef CGAL::Homogeneous<testnt> TestR;
typedef TestR::RT inputt;
#endif

#if TESTR == 5

typedef leda_rational testnt;
typedef CGAL::Cartesian<testnt> TestR;
typedef int inputt;
#endif


#if TESTR == 6

typedef leda_rational testnt;
typedef CGAL::Homogeneous<testnt> TestR;
typedef int inputt;
#endif

#if TESTR == 7

typedef CGAL::TestfieldC testnt;
typedef CGAL::Cartesian<testnt> TestR;
typedef TestR::RT inputt;

#endif


#if TESTR == 8

typedef CGAL::TestrepH testnt;
typedef CGAL::Homogeneous<testnt> TestR;
typedef TestR::RT inputt;

#endif

#if BBOX == 1
typedef CGAL::Pm_segment_traits_2<TestR>        Traits;
#else
typedef CGAL::Pm_straight_traits_2<TestR>       Traits;
#endif
