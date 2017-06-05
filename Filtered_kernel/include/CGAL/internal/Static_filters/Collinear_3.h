// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {


#include <iostream>


template < typename K_base >
class Collinear_3
  : public K_base::Collinear_3
{
  typedef typename K_base::Point_3          Point_3;
  typedef typename K_base::Vector_3         Vector_3;
  typedef typename K_base::Sphere_3         Sphere_3;
  typedef typename K_base::Tetrahedron_3    Tetrahedron_3;
  typedef typename K_base::Collinear_3    Base;

public:
 typedef typename Base::result_type  result_type;
  result_type 
  operator()(const Point_3 &p, const Point_3 &q,
	     const Point_3 &r, const Point_3 &s) const
  {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Collinear_3", tmp);

      Get_approx<Point_3> get_approx; // Identity functor for all points
                                      // but lazy points.

      double px, py, pz, qx, qy, qz, rx, ry, rz;

      if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
          fit_in_double(get_approx(p).z(), pz) &&
          fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
          fit_in_double(get_approx(q).z(), qz) &&
          fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
          fit_in_double(get_approx(r).z(), rz))
        {

        }
      return Base::operator()(p, q, r);
  }

};

    } } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COLLINEAR_3_H
