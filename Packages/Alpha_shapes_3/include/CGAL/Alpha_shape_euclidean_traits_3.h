// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/basic.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/smallest_radius_3.h>

CGAL_BEGIN_NAMESPACE

//------------------ Function Objects ---------------------------------

template < class return_type, class T >
class Compute_squared_radius_circumsphere_3
{
public:

  typedef return_type result_type;

  result_type operator()(const T& p, const T& q, const T& r, const T& s) const
    { 
      return CGAL::squared_radius(p, q, r, s);
    }

  result_type operator()(const T& p, const T& q, const T& r) const
    { 
      return CGAL::squared_radius(p, q, r);
    }

  result_type operator()(const T& p, const T& q) const
    { 
      return CGAL::squared_radius_smallest_circumsphere(p, q);
    }
};

//------------------ Traits class -------------------------------------

template <class R>
class Alpha_shape_euclidean_traits_3 : public R
{
public:
 
  typedef R Rep;
  typedef typename R::FT Coord_type;
  typedef typename R::Point_3 Point_3;
  typedef Point_3  Point;
  typedef typename R::Segment_3 Segment_3;

  typedef CGAL::Compute_squared_radius_circumsphere_3<Coord_type, Point_3> 
  Compute_squared_radius_circumsphere_3;
  typedef typename R::Side_of_bounded_sphere_3 Side_of_bounded_sphere_3;


  Compute_squared_radius_circumsphere_3 compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_circumsphere_3();
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
