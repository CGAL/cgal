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
// file          : include/CGAL/Alpha_shape_euclidean_traits_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef ALPHA_SHAPE_TRAITS_H
#define ALPHA_SHAPE_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
//#include <CGAL/Homogeneous.h>
//#include <CGAL/Integer.h>
//#include <CGAL/Rational.h>
//#include <CGAL/Fixed.h>
//#include <CGAL/Real.h>

#include <CGAL/squared_distance_3.h>   // to avoid a g++ problem
#include <CGAL/Point_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>

#include <CGAL/Smallest_radius_3.h>
#include <CGAL/Predicates_on_points_3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template<class FT>
VectorC3<FT> operator^(const VectorC3<FT> &v,const VectorC3<FT> &w)
  {
    return VectorC3<FT>( v.y() * w.z() - v.z() * w.y() ,
			 v.z() * w.x() - v.x() * w.z() ,
			 v.x() * w.y() - v.y() * w.x() );
  }

//------------------ Traits class -------------------------------------

template<class R>
class Alpha_shape_euclidean_traits_3 : public Triangulation_geom_traits_3<R>
{
public: 
  typedef R::FT Coord_type;

 //---------------------------------------------------------------------
  Coord_type squared_radius_smallest_circumsphere(const Point_3<R> &p,
						  const Point_3<R> &q,
						  const Point_3<R> &r,
						  const Point_3<R> &s) const
  {
    return squared_radius_smallest_circumsphere(p, q, r, s);
  }

  //---------------------------------------------------------------------

  Coord_type squared_radius_smallest_circumsphere(const Point_3<R> &p,
						  const Point_3<R> &q,
						  const Point_3<R> &r) const
  {
    return squared_radius_smallest_circumsphere(p, q, r);
  }
  
  //------------------------------------------------------------------

  Coord_type squared_radius_smallest_circumsphere(const Point_3<R> &p,
						  const Point_3<R> &q) const
  {
    return squared_radius_smallest_circumsphere(p, q);
  }

  //------------------------------------------------------------------

  Bounded_side side_of_bounded_sphere(const Point_3<R> &p,
					   const Point_3<R> &q,
					   const Point_3<R> &test) const
  {
    return side_of_bounded_sphere(p, q, test);
  }

  //------------------------------------------------------------------

  Bounded_side side_of_bounded_sphere(const Point_3<R> &p,
					   const Point_3<R> &q,
					   const Point_3<R> &r,
					   const Point_3<R> &test) const
  {
    return side_of_bounded_sphere(p, q, r, test);
  }
  //------------------------------------------------------------------
  //EXEMPLE D'ANDREAS
  Squared_radius_smallest_circumsphere
  squared_radius_smallest_circumsphere_object() const {
    return Squared_radius_smallest_circumsphere();
  }

};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif
