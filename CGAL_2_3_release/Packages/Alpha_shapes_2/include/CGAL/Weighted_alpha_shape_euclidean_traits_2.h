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
// file          : include/CGAL/Weighted_alpha_shape_euclidean_traits_2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H

#include <CGAL/basic.h>

#include <CGAL/squared_distance_2.h>

#include <CGAL/in_smallest_orthogonalcircleC2.h>
#include <CGAL/squared_radius_smallest_orthogonalcircleC2.h>

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

//------------------ Function Objects----------------------------------

template < class return_type, class T >
class Compute_squared_radius_orthogonalcircle_2
{
public:
  typedef return_type result_type;
  
  result_type operator()(const T& p, const T& q, const T& r)
    {
      return
	std::max
	(return_type(0), CGAL::squared_radius_orthogonalcircle(p, q, r));
    }

  result_type operator()(const T& p, const T& q)
    {
      return
	std::max
	(return_type(0), CGAL::squared_radius_smallest_orthogonalcircle(p, q));
    }
};

//-------------------------------------------------------------------

template < class T >
class Side_of_bounded_orthogonalcircle_2
{
public:
  typedef Bounded_side result_type;
  
   result_type operator()(const T& p, const T& q, const T& t)
    {  
      return
	CGAL::in_smallest_orthogonalcircle(p, q, t);
    }
};

//------------------ Traits class -------------------------------------

template< class R >
class Weighted_alpha_shape_euclidean_traits_2 : public
Regular_triangulation_euclidean_traits_2<R, typename R::FT> 
{

public: 
  
  typedef typename R::FT Coord_type;
  typedef typename 
   Regular_triangulation_euclidean_traits_2<R, typename R::FT>::Weighted_point 
     Point;

  typedef Compute_squared_radius_orthogonalcircle_2<Coord_type, Point>
  Compute_squared_radius_orthogonalcircle_2;
  typedef Side_of_bounded_orthogonalcircle_2<Point>
  Side_of_bounded_orthogonalcircle_2;
  
  //------------------------------------------------------------------

  Compute_squared_radius_orthogonalcircle_2
  compute_squared_radius_2_object() const
    {
      return Compute_squared_radius_orthogonalcircle_2();
    }

  //------------------------------------------------------------------

  Side_of_bounded_orthogonalcircle_2 side_of_bounded_circle_2_object() const
    {
      return Side_of_bounded_orthogonalcircle_2();
    }
};
  
//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif
