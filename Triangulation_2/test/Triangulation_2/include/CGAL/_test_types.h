#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <CGAL/_test_short_names_2.h>

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>




// #ifdef CGAL_USE_CLN
// #include <CGAL/CLN/cl_integer.h>
// #include <cl_io.h>
// #include <cl_integer_io.h>
// typedef cl_I Rtype;
// #else
// #ifdef CGAL_USE_GMP
// #include <CGAL/Gmpz.h>
// typedef CGAL::Gmpz Rtype;
// #else
// #ifdef CGAL_USE_LEDA
// #include <CGAL/leda_integer.h>
// typedef leda_integer Rtype;
// #else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float Rtype;
//#endif // CGAL_USE_LEDA
//#endif // CGAL_USE_GMP
//#endif // CGAL_USE_CLN

#include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float>  Ftype;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>

typedef CGAL::Simple_cartesian<Ftype>         Test_rep_cartesian;
typedef CGAL::Homogeneous<Rtype>              Test_rep_homogeneous;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       K;
struct TestK : public K {}; 

typedef CGAL::Quotient<Ftype>                       Exact_type;
typedef CGAL::Simple_cartesian<Exact_type>          Exact_kernel;
struct EK : public Exact_kernel {};

#endif
