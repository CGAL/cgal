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

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/constructions/constructions_on_weighted_points_cartesian_3.h>
#include <CGAL/predicates/predicates_on_weighted_points_cartesian_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

CGAL_BEGIN_NAMESPACE

//------------------ Function Objects----------------------------------

template < class K >
class Compute_squared_radius_orthogonal_sphere_3
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
    FT res = squared_radius_orthogonal_sphereC3(px, py, pz, pw,
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
    
    FT res = squared_radius_smallest_orthogonal_sphereC3(px, py, pz, pw,
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

    FT res = squared_radius_smallest_orthogonal_sphereC3(px, py, pz, pw,
						    qx, qy, qz, qw);
    return max (FT(0), res);
    }
};


   
//------------------ Traits class -------------------------------------

template <class R>
class Weighted_alpha_shape_euclidean_traits_3 : public 
Regular_triangulation_euclidean_traits_3<R>
{
public:
 
  typedef Weighted_alpha_shape_euclidean_traits_3<R> Self;
  typedef Regular_triangulation_euclidean_traits_3<R> Base;
  typedef R Rep;
  typedef typename R::FT Coord_type;

  typedef typename Base::Bare_point      Bare_Point;
  typedef typename Base::Weighted_point  Weighted_point;
  typedef Weighted_point Point_3;
  typedef Weighted_point Point;

  typedef CGAL::Compute_squared_radius_orthogonal_sphere_3<Self> 
    Compute_squared_radius_orthogonal_sphere_3;
  typedef CGAL::Side_of_bounded_orthogonal_sphere_3<R> 
    Side_of_bounded_orthogonal_sphere_3;

  //---------------------------------------------------------------------

  Compute_squared_radius_orthogonal_sphere_3 
  compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_orthogonal_sphere_3();
    }
  //---------------------------------------------------------------------

  Side_of_bounded_orthogonal_sphere_3 
  side_of_bounded_sphere_3_object() const
    {
      return Side_of_bounded_orthogonal_sphere_3();
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 
