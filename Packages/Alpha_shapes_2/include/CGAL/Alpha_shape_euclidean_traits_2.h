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
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H
#define CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
//#include <CGAL/Homogeneous.h>
//#include <CGAL/Integer.h>
//#include <CGAL/Rational.h>
//#include <CGAL/Fixed.h>

#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/smallest_radius_2.h>
#include <CGAL/side_of_smallest_circle_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

//------------------ Traits class -------------------------------------

template<class R>
class Alpha_shape_euclidean_traits_2 : public
Triangulation_euclidean_traits_2<R> 
{
public: 
  typedef typename R::FT Coord_type;

  //---------------------------------------------------------------------

  Coord_type squared_radius_smallest_circumcircle(const Point &p,
						  const Point &q,
						  const Point &r)
    const 
    {
      // the computation of the squared radius takes 17 multiplications
      // and 12 additions

      return CGAL::squared_radius_smallest_circumcircle(p, q, r);
    }
  
  //------------------------------------------------------------------

  Coord_type squared_radius_smallest_circumcircle(const Point &p,
						  const Point &q)
    const 
    {
      // the computation of the squared radius takes 17 multiplications
      // and 12 additions

      return CGAL::squared_radius_smallest_circumcircle(p, q);

    }

  //------------------------------------------------------------------

  Bounded_side side_of_bounded_circle(const Point &p,
				      const Point &q,
				      const Point &test) const 
    {
      return CGAL::side_of_bounded_circle(p, q, test);
    }

  //------------------------------------------------------------------

};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif
