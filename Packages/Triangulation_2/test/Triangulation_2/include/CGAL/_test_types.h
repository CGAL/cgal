#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <iostream.h>
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <assert.h>
#ifndef CGAL_GMPZ_H
#include <CGAL/Gmpz.h>
#endif // CGAL_GMPZ
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H

#include <CGAL/_test_short_names_2.h>

typedef CGAL_Cartesian<CGAL_Gmpz> Test_rep_cartesian;
typedef CGAL_Homogeneous<CGAL_Gmpz> Test_rep_homogeneous;

//the following define shorter names to make the (g++/egcs) linker work 
//#define CGAL__Triangulation_test_traits CGAL_Ttt

 class CGAL__Triangulation_test_traits;

#endif
