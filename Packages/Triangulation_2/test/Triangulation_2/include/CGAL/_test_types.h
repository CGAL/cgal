#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <CGAL/_test_short_names_2.h>

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <iostream>
#include <cassert>
// #ifndef CGAL_GMPZ_H
// #include <CGAL/Gmpz.h>
// #endif // CGAL_GMPZ_H
// #include <CGAL/Quotient.h>
//#include<CGAL/Arithmetic_filter.h>
#include <CGAL/leda_integer.h>


#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H


// typedef CGAL::Quotient<CGAL::Gmpz> CT;
// typedef CGAL::Quotient<CGAL::Gmpz> ET;
// typedef CGAL::Filtered_exact<CT,ET> coord_type;
// typedef CGAL::Quotient<CGAL::Gmpz> coord_type;
// typedef CGAL::Cartesian<coord_type> Test_rep_cartesian;
// typedef CGAL::Homogeneous<CGAL::Gmpz> Test_rep_homogeneous;
typedef leda_integer                   Rtype;
typedef CGAL::Quotient<leda_integer>   Ftype;
typedef CGAL::Cartesian<Ftype>         Test_rep_cartesian;
typedef CGAL::Homogeneous<Rtype>       Test_rep_homogeneous;

#endif
