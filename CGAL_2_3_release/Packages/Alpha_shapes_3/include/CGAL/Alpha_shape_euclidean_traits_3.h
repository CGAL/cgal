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

#ifndef CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

//#include <CGAL/basic.h>


#include <CGAL/squared_distance_3.h>   // to avoid a g++ problem
#include <CGAL/Point_3.h>
#include <CGAL/predicates_on_points_3.h>

#include <CGAL/Triangulation_geom_traits_3.h>

#include <CGAL/smallest_radius_3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

//------------------ Function Objects ---------------------------------

template < class return_type, class T >
class Compute_squared_radius_circumsphere_3
{
public:

  typedef return_type result_type;

  result_type operator()(const T& p, const T& q, const T& r, const T& s)
    { 
      return CGAL::squared_radius(p, q, r, s);
    }

  result_type operator()(const T& p, const T& q, const T& r)
    { 
      return CGAL::squared_radius(p, q, r);
    }

  result_type operator()(const T& p, const T& q)
    { 
      return CGAL::squared_radius_smallest_circumsphere(p, q);
    }
};

//------------------ Traits class -------------------------------------

template <class R>
class Alpha_shape_euclidean_traits_3 : public Triangulation_geom_traits_3<R>
{
public:
 
  typedef R Rep;
  typedef typename R::FT Coord_type;
  typedef typename Triangulation_geom_traits_3<R>::Point_3 Point_3;
  typedef Point_3  Point;
  typedef typename Triangulation_geom_traits_3<R>::Segment_3 Segment_3;

  typedef Compute_squared_radius_circumsphere_3<Coord_type, Point_3> 
  Compute_squared_radius_circumsphere_3;
  typedef typename R::Side_of_bounded_sphere_3 Side_of_bounded_sphere_3;

  //---------------------------------------------------------------------

  Compute_squared_radius_circumsphere_3 compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_circumsphere_3();
    }
  //---------------------------------------------------------------------
  
  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object() const
    {
      return Side_of_bounded_sphere_3();
    }
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
