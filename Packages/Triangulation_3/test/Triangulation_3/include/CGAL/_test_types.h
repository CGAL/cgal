#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <iostream>
#include <cassert>

// #ifndef CGAL_GMPZ_H
// #include <CGAL/Gmpz.h>
// #endif // CGAL_GMPZ_H

// #ifdef CGAL_USE_CLN
// #include <CGAL/CLN/cl_integer.h>
// typedef cl_I my_NT;
// #else
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer my_NT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz my_NT;
#else
#include <CGAL/double.h>
typedef double my_NT;
#endif
#endif
// #endif

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H

// #include <CGAL/_test_short_names_2.h>

typedef CGAL::Cartesian<my_NT> Test_rep_cartesian;
// #include <CGAL/leda_real.h>
// typedef CGAL::Homogeneous<my_NT,leda_real> Test_rep_homogeneous;
typedef CGAL::Homogeneous<my_NT> Test_rep_homogeneous;

//the following define shorter names to make the (g++/egcs) linker work 
// #define CGAL__Triangulation_test_traits_3 CGAL::Ttt3

// class CGAL__Triangulation_test_traits;
class _Triangulation_test_traits;

#endif
