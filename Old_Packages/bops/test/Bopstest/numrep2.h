#if TestNumCode == 1
typedef float TestNum;
#endif

#if TestNumCode == 2
typedef double TestNum;
#endif

#if TestNumCode == 3
typedef long TestNum;
#endif

#if TestNumCode == 4
#include <CGAL/leda_integer.h>
typedef leda_integer TestNum;
#endif

#if TestNumCode == 5
#include <CGAL/leda_rational.h>
typedef leda_rational TestNum;
#endif

#if TestNumCode == 6
typedef int TestNum;
#endif

#if TestNumCode == 7
typedef CGAL::Gmpz TestNum;
#endif


#if TestRepCode == 2
typedef CGAL::Cartesian<TestNum> TestRep;
#else
typedef CGAL::Homogeneous<TestNum> TestRep;
#endif
