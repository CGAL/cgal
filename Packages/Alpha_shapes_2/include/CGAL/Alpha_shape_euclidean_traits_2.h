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
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H
#define CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H

#include <CGAL/basic.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>


//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

//------------------ Function Objects ---------------------------------

template < class K>
class Compute_squared_radius_2
{
public:
  typedef typename K::FT result_type;
  typedef typename K::Point_2 Point;

  result_type operator()(const Point& p, const Point& q, const Point& r) const
    {
      return CGAL::squared_radius(p, q, r);
    }

  result_type operator()(const Point& p, const Point& q) const
    {
      typename K::Vector_2 v(p-q);
      return (v*v)/4;
    }
};

//------------------ Traits class -------------------------------------

template < class R >
class Alpha_shape_euclidean_traits_2 : public
Triangulation_euclidean_traits_2<R> 
{
public: 
  typedef typename R::FT Coord_type;
  typedef typename R::Point_2 Point;

  typedef CGAL::Compute_squared_radius_2<R> Compute_squared_radius_2;
  typedef typename R::Side_of_bounded_circle_2 Side_of_bounded_circle_2;
  
  //------------------------------------------------------------------

  Compute_squared_radius_2
  compute_squared_radius_2_object() const
    {
      return Compute_squared_radius_2();
    }

  //------------------------------------------------------------------

  Side_of_bounded_circle_2 side_of_bounded_circle_2_object() const
    {
      return Side_of_bounded_circle_2();
    }
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif
