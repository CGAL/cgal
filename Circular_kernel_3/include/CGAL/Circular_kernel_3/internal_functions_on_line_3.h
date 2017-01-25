// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {
  namespace SphericalFunctors {

    template < class SK >
    typename SK::Line_3
    construct_line_3(const typename SK::Polynomials_for_line_3 &eq)
    {
      typedef typename SK::Line_3 Line_3;
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Vector_3 Vector_3;
      return Line_3(Point_3(eq.b1(),eq.b2(),eq.b3()),
                    Vector_3(eq.a1(),eq.a2(),eq.a3()));
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_LINE_3_H
