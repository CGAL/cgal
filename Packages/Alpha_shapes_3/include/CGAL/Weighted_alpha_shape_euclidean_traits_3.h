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
// file          : include/CGAL/Weighted_alpha_shape_euclidean_traits_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/constructions/squared_radius_smallest_orthogonalsphere_ftC3.h>
#include <CGAL/predicates/in_smallest_orthogonalsphere_ftC3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

CGAL_BEGIN_NAMESPACE

//------------------ Function Objects----------------------------------

template < class K >
class Compute_squared_radius_orthogonalsphere_3
{
public:

  typedef typename K::Point Point;
  typedef typename K::FT FT;
  typedef typename K::FT return_type;
  typedef Arity_tag< 4 >   Arity;

  return_type operator()(const Point& p, const Point& q, 
			 const Point& r, const Point& s) const
  {   
    FT px(p.point().x());
    FT py(p.point().y());
    FT pz(p.point().z());
    FT pw(p.weight());
    FT qx(q.point().x());
    FT qy(q.point().y());
    FT qz(q.point().z());
    FT qw(q.weight());
    FT rx(r.point().x());
    FT ry(r.point().y()); 
    FT rz(r.point().z());
    FT rw(r.weight()); 
    FT sx(s.point().x());
    FT sy(s.point().y());
    FT sz(s.point().z());
    FT sw(s.weight());
    FT res = squared_radius_orthogonalsphereC3(px, py, pz, pw,
					       qx, qy, qz, qw,
					       rx, ry, rz, rw,
					       sx, sy, sz, sw);
      return max (FT(0), res);
    }

  return_type operator()(const Point& p, const Point& q, const Point& r) const
  {  
    FT px(p.point().x());
    FT py(p.point().y());
    FT pz(p.point().z());
    FT pw(p.weight());
    FT qx(q.point().x());
    FT qy(q.point().y());
    FT qz(q.point().z());
    FT qw(q.weight());
    FT rx(r.point().x());
    FT ry(r.point().y()); 
    FT rz(r.point().z());
    FT rw(r.weight()); 
    
    FT res = squared_radius_smallest_orthogonalsphereC3(px, py, pz, pw,
							qx, qy, qz, qw,
							rx, ry, rz, rw); 
    return max (FT(0), res );
  }

  return_type operator()(const Point& p, const Point& q) const
  {   
    FT px(p.point().x());
    FT py(p.point().y());
    FT pz(p.point().z());
    FT pw(p.weight());
    FT qx(q.point().x());
    FT qy(q.point().y());
    FT qz(q.point().z());
    FT qw(q.weight());

    res = squared_radius_smallest_orthogonalsphereC3(px, py, pz, pw,
						    qx, qy, qz, qw);
    return max (FT(0), res);
    }
};

//-------------------------------------------------------------------

template < class K >
class Side_of_bounded_orthogonalsphere_3
{
public:
  typedef typename K::Point Point;
  typedef typename K::FT FT;
  typedef Bounded_side result_type;
  typedef Arity_tag< 4 >   Arity;

  result_type
  operator()(const Point& p, const Point& q, const Point& r,
	     const Point& t) const
  {  
    FT px(p.point().x());
    FT py(p.point().y());
    FT pz(p.point().z());
    FT pw(p.weight());
    FT qx(q.point().x());
    FT qy(q.point().y());
    FT qz(q.point().z());
    FT qw(q.weight());
    FT rx(r.point().x());
    FT ry(r.point().y());
    FT rz(r.point().z());
    FT rw(r.weight());
    FT tx(t.point().x());
    FT ty(t.point().y());
    FT tz(t.point().z());
    FT tw(t.weight());
    
    return in_smallest_orthogonalsphereC3(px, py, pz, pw,
					  qx, qy, qz, qw,
					  rx, ry, rz, rw,
					  tx, ty, tz, tw);
  }
};
  
//------------------ Traits class -------------------------------------

template <class R>
class Weighted_alpha_shape_euclidean_traits_3 : public 
Regular_triangulation_euclidean_traits_3<R>
{
public:
 
  typedef Weighted_alpha_shape_euclidean_traits_3<R> Self;
  typedef R Rep;
  typedef typename R::FT Coord_type;

  typedef typename 
    Regular_triangulation_euclidean_traits_3<R>::Bare_point Bare_Point;
  typedef typename 
    Regular_triangulation_euclidean_traits_3<R>::Weighted_point Weighted_point;
  typedef Weighted_point Point_3;
  typedef Weighted_point Point;

  typedef CGAL::Compute_squared_radius_orthogonalsphere_3<Self> 
    Compute_squared_radius_orthogonalsphere_3;
  typedef CGAL::Side_of_bounded_orthogonalsphere_3<Self> 
    Side_of_bounded_orthogonalsphere_3;

  //---------------------------------------------------------------------

  Compute_squared_radius_orthogonalsphere_3 
  compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_orthogonalsphere_3();
    }
  //---------------------------------------------------------------------

  Side_of_bounded_orthogonalsphere_3 
  side_of_bounded_sphere_3_object() const
    {
      return Side_of_bounded_orthogonalsphere_3();
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 
