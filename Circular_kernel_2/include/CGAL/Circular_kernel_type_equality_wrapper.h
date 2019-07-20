// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)


#ifndef CGAL_CIRCULAR_KERNEL_TYPE_EQUALITY_WRAPPER_H
#define CGAL_CIRCULAR_KERNEL_TYPE_EQUALITY_WRAPPER_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/user_classes.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Circular_arc_point_2.h>
#include <CGAL/Line_arc_2.h>
#include <CGAL/Root_of_traits.h>


namespace CGAL {

// This is a kernel wrapper which provides the type equality between
// Kernel::Point_2 and CGAL::Point_2<Kernel>, by deriving from
// K_base::Point_2 (and similar for the other types).

template < typename K_base, typename Kernel >
struct Circular_kernel_type_equality_wrapper
  : public Type_equality_wrapper<K_base, Kernel>
{
    typedef K_base                                  Kernel_base;
    typedef CGAL::Circular_arc_2<Kernel>               Circular_arc_2;     
    typedef CGAL::Line_arc_2<Kernel>                   Line_arc_2;
    typedef CGAL::Circular_arc_point_2<Kernel>         Circular_arc_point_2;

    //Something has to be done with these 3, maybe a lazy Algebraic kernel?
	   
    //typedef Polynomial_for_circles_2_2<Kernel>   Polynomial_for_circles_2_2;
    //typedef Polynomial_1_2<Kernel>               Polynomial_1_2;
};

} //namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_TYPE_EQUALITY_WRAPPER_H
