#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <CGAL/_test_short_names_2.h>

#include <iostream>
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <cassert>
#ifndef CGAL_GMPZ_H
#include <CGAL/Gmpz.h>
#endif // CGAL::GMPZ
#include <CGAL/Quotient.h>
//#include<CGAL/Arithmetic_filter.h>

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H


typedef CGAL::Quotient<CGAL::Gmpz> CT;
typedef CGAL::Quotient<CGAL::Gmpz> ET;
//typedef CGAL::Filtered_exact<CT,ET> coord_type;
typedef CGAL::Quotient<CGAL::Gmpz> coord_type;
typedef CGAL::Cartesian<coord_type> Test_rep_cartesian;
//typedef CGAL::Cartesian<CGAL::Gmpz>  Test_rep_cartesian;
typedef CGAL::Homogeneous<CGAL::Gmpz> Test_rep_homogeneous;



#endif
