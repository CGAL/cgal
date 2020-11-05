#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <iostream>
#include <cassert>

#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float Rtype;

#include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float>  Ftype;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>

typedef CGAL::Simple_cartesian<Ftype>         Test_rep_cartesian;
typedef CGAL::Homogeneous<Rtype>              Test_rep_homogeneous;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1> TestK;

typedef CGAL::Quotient<Ftype>                       Exact_type;
typedef CGAL::Simple_cartesian<Exact_type>          EK;

#endif
