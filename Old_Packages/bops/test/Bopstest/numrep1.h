// TestNumCode should be either
// 1: float
// 2: double
// 3: long
// 4: Integer
// 5: Rational
// 6: int

// TestRepCode should be either
// 1: Homogeneous
// 2: Cartesian

#include <CGAL/number_type_tags.h>

#ifndef TestNumCode
#   define TestNumCode 7
#endif

#ifndef TestRepCode
#   define TestRepCode 1
#endif

#if TestRepCode == 2
#include <CGAL/Cartesian.h>
#else
#include <CGAL/Homogeneous.h>
//#include <CGAL/Cartesian.h> // because of bug in bops code.
#endif


#if TestNumCode == 4
#include <CGAL/leda_integer.h>
#endif

#if TestNumCode == 5
#include <CGAL/leda_rational.h>
#endif

#if TestNumCode == 6
#endif

#if TestNumCode == 7
#include <CGAL/Gmpz.h>
#endif

