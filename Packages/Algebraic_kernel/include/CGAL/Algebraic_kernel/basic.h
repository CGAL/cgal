// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_BASIC_H
#define CGAL_ALGEBRAIC_KERNEL_BASIC_H

// This file currently gathers function declarations.
// This is needed by the way the current NT interface works.

#include <utility>
#include <CGAL/config.h>

CGAL_BEGIN_NAMESPACE

class MP_Float;
template < typename T > class Root_of_2;
template < typename T > class Lazy_exact_nt;

template < typename T >
std::pair<double,double> to_interval(const Root_of_2<T>&);

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_BASIC_H
