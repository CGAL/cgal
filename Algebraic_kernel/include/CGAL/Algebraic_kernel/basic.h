// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion

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

template < typename RT >
Root_of_2<RT> square(const Root_of_2<RT> &);

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_BASIC_H
