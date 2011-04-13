#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <CGAL/_test_short_names_2.h>

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <iostream>
#include <cassert>

// #ifdef CGAL_USE_CLN
// #include <CGAL/CLN/cl_integer.h>
// #include <cl_io.h>
// #include <cl_integer_io.h>
// typedef cl_I Rtype;
// #else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz Rtype;
#else
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer Rtype;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float Rtype;
#endif // CGAL_USE_LEDA
#endif // CGAL_USE_GMP
//#endif // CGAL_USE_CLN

#include <CGAL/Quotient.h>
typedef CGAL::Quotient<Rtype>   Ftype;

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>

typedef CGAL::Cartesian<Ftype>         Test_rep_cartesian;
typedef CGAL::Homogeneous<Rtype>       Test_rep_homogeneous;

// typedef CGAL::Quotient<CGAL::Gmpz> CT;
// typedef CGAL::Quotient<CGAL::Gmpz> ET;
// typedef CGAL::Filtered_exact<CT,ET> coord_type;
// typedef CGAL::Quotient<CGAL::Gmpz> coord_type;
// typedef CGAL::Cartesian<coord_type> Test_rep_cartesian;
// typedef CGAL::Homogeneous<CGAL::Gmpz> Test_rep_homogeneous;

#endif
