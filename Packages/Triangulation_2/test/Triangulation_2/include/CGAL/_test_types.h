#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <iostream>
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <cassert>
#ifndef CGAL_GMPZ_H
#include <CGAL/Gmpz.h>
#endif // CGAL::GMPZ
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H

#include <CGAL/_test_short_names_2.h>

typedef CGAL::Cartesian<CGAL::Gmpz> Test_rep_cartesian;
typedef CGAL::Homogeneous<CGAL::Gmpz> Test_rep_homogeneous;

//the following define shorter names to make the (g++/egcs) linker work 
//#define CGAL::_Triangulation_test_traits CGAL::Ttt

 class _Triangulation_test_traits;

#endif
