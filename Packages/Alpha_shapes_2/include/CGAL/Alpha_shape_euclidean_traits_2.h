// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/Alpha_shape_euclidean_traits_2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

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
