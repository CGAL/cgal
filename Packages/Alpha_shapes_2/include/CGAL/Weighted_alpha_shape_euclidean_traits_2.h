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
#include <CGAL/Cartesian.h>

//#include <CGAL/Homogeneous.h>
//#include <CGAL/Integer.h>
//#include <CGAL/Rational.h>
//#include <CGAL/Fixed.h>

#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/in_smallest_orthogonalcircleC2.h>
#include <CGAL/squared_radius_smallest_orthogonalcircleC2.h>

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

//------------------ Traits class -------------------------------------

template< class R >
class Weighted_alpha_shape_euclidean_traits_2 : public
Regular_triangulation_euclidean_traits_2<R, typename R::FT> 
{

public: 
  
  typedef typename R::FT Coord_type;
  typedef typename 
  Regular_triangulation_euclidean_traits_2<R, typename R::FT>::Weighted_point Point;

  //---------------------------------------------------------------------

  Coord_type squared_radius(const Point &p,
			    const Point &q,
			    const Point &r) const 
    {
  
      return
	std::max
	(Coord_type(0), CGAL::squared_radius_orthogonalcircle(p, q, r));
    }




  Coord_type squared_radius(const Point &p,
			    const Point &q) const 
    {

      return
	std::max
	(Coord_type(0), CGAL::squared_radius_smallest_orthogonalcircle(p, q));
    }

  Bounded_side side_of_circle(const Point &p,
			      const Point &q,
			      const Point &t) const 
    {
  
      return
	CGAL::in_smallest_orthogonalcircle(p, q, t);
    }

};
  
//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif
