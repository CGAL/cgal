// Copyright (c) 1998  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri


#ifndef CGAL_DISTANCE_3_0_H
#define CGAL_DISTANCE_3_0_H

#include <CGAL/kernel_assertions.h>
#include <CGAL/enum.h>
#include <CGAL/wmult.h>

#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>


namespace CGAL {


namespace internal {

template <class K>
bool is_null(const typename K::Vector_3 &v, const K&)
{
    typedef typename K::RT RT;
    return v.hx()==RT(0) && v.hy()==RT(0) && v.hz()==RT(0);
}



template <class K>
typename K::RT
wdot(const typename K::Vector_3 &u,
     const typename K::Vector_3 &v, 
     const K&)
{
    return  (u.hx()*v.hx() + u.hy()*v.hy() + u.hz()*v.hz());
}

template <class K>
typename K::RT
wdot_tag(const typename K::Point_3 &p,
	 const typename K::Point_3 &q,
	 const typename K::Point_3 &r, 
	 const K&,
	 const Cartesian_tag&)
{
  return  ( (p.x() - q.x()) * (r.x()- q.x())
	    + (p.y()- q.y())* (r.y()- q.y())
	    + (p.z()- q.z()) * (r.z() - q.z()) );
}

template <class K>
typename K::RT
wdot_tag(const typename K::Point_3 &p,
	 const typename K::Point_3 &q,
	 const typename K::Point_3 &r, 
	 const K&,
	 const Homogeneous_tag&)
{
    return  ( (p.hx() * q.hw() - q.hx() * p.hw())
	      * (r.hx() * q.hw() - q.hx() * r.hw())
	      + (p.hy() * q.hw() - q.hy() * p.hw())
	      * (r.hy() * q.hw() - q.hy() * r.hw())
	      + (p.hz() * q.hw() - q.hz() * p.hw())
	      * (r.hz() * q.hw() - q.hz() * r.hw()));
}

template <class K>
inline
typename K::RT
wdot(const typename K::Point_3 &p,
     const typename K::Point_3 &q,
     const typename K::Point_3 &r, 
     const K& k)
{ 
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return wdot_tag(p, q, r, k, tag);
}




template <class K>
typename K::Vector_3
wcross(const typename K::Vector_3 &u,
       const typename K::Vector_3 &v, 
       const K&)
{
  typedef typename K::Vector_3 Vector_3;
    return Vector_3(
        u.hy()*v.hz() - u.hz()*v.hy(),
        u.hz()*v.hx() - u.hx()*v.hz(),
        u.hx()*v.hy() - u.hy()*v.hx());
}


template <class K>
inline 
bool 
is_acute_angle(const typename K::Vector_3 &u,
	       const typename K::Vector_3 &v, 
	       const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) > RT(0) ;
}

template <class K>
inline 
bool 
is_straight_angle(const typename K::Vector_3 &u,
		  const typename K::Vector_3 &v, 
		  const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) == RT(0) ;
}

template <class K>
inline 
bool 
is_obtuse_angle(const typename K::Vector_3 &u,
		const typename K::Vector_3 &v, 
		const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) < RT(0) ;
}

template <class K>
inline 
bool 
is_acute_angle(const typename K::Point_3 &p,
	       const typename K::Point_3 &q,
	       const typename K::Point_3 &r, 
	       const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) > RT(0) ;
}

template <class K>
inline
bool 
is_straight_angle(const typename K::Point_3 &p,
		  const typename K::Point_3 &q,
		  const typename K::Point_3 &r, 
		  const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) == RT(0) ;
}

template <class K>
inline 
bool 
is_obtuse_angle(const typename K::Point_3 &p,
		const typename K::Point_3 &q, 
		const typename K::Point_3 &r, 
		  const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) < RT(0) ;
}
template <class K>
inline 
typename K::FT
squared_distance(const typename K::Point_3 & pt1,
		 const typename K::Point_3 & pt2, 
		 const K& k)
{
  return k.compute_squared_distance_3_object()(pt1, pt2);
}


template <class K>
typename K::FT
squared_distance_to_plane(const typename K::Vector_3 & normal,
			  const typename K::Vector_3 & diff, 
			  const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  RT dot, squared_length;
  dot = wdot(normal, diff, k);
  squared_length = wdot(normal, normal, k);
  return FT(dot*dot) /
    FT(wmult((K*)0, squared_length, diff.hw(), diff.hw()));
}


template <class K>
typename K::FT
squared_distance_to_line(const typename K::Vector_3 & dir,
			 const typename K::Vector_3 & diff, 
			 const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  Vector_3 wcr = wcross(dir, diff, k);
  return FT(wcr*wcr)/FT(wmult(
			      (K*)0, RT(wdot(dir, dir, k)), diff.hw(), diff.hw()));
}


template <class K>
inline
bool
same_direction_tag(const typename K::Vector_3 &u,
		   const typename K::Vector_3 &v,
		   const K&,
		   const Cartesian_tag&)
{ 
  typedef typename K::FT FT;
  const FT& ux = u.x();
  const FT& uy = u.y();
  const FT& uz = u.z();
  if (CGAL_NTS abs(ux) > CGAL_NTS abs(uy)) {
    if (CGAL_NTS abs(ux) > CGAL_NTS abs(uz)) {
      return CGAL_NTS sign(ux) == CGAL_NTS sign(v.x());
    } else {
      return CGAL_NTS sign(uz) == CGAL_NTS sign(v.z());
    }
  } else {
    if (CGAL_NTS abs(uy) > CGAL_NTS abs(uz)) {
      return CGAL_NTS sign(uy) == CGAL_NTS sign(v.y());
    } else {
      return CGAL_NTS sign(uz) == CGAL_NTS sign(v.z());
    }
  } 
}


template <class K>
inline
bool
same_direction_tag(const typename K::Vector_3 &u,
		   const typename K::Vector_3 &v,
		   const K&,
		   const Homogeneous_tag&)
{   
  typedef typename K::RT RT;
  const RT& uhx = u.hx();
  const RT& uhy = u.hy();
  const RT& uhz = u.hz();
  if (CGAL_NTS abs(uhx) > CGAL_NTS abs(uhy)) {
    if (CGAL_NTS abs(uhx) > CGAL_NTS abs(uhz)) {
      return CGAL_NTS sign(uhx) == CGAL_NTS sign(v.hx());
    } else {
      return CGAL_NTS sign(uhz) == CGAL_NTS sign(v.hz());
    }
  } else {
    if (CGAL_NTS abs(uhy) > CGAL_NTS abs(uhz)) {
      return CGAL_NTS sign(uhy) == CGAL_NTS sign(v.hy());
    } else {
      return CGAL_NTS sign(uhz) == CGAL_NTS sign(v.hz());
    }
  }
}


template <class K>
inline
bool
same_direction(const typename K::Vector_3 &u,
	       const typename K::Vector_3 &v,
	       const K& k)
{  
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return same_direction_tag(u, v, k, tag);
}


} // namespace internal

template <class K>
inline 
typename K::FT
squared_distance(const Point_3<K> & pt1,
		 const Point_3<K> & pt2)
{
  return internal::squared_distance(pt1,pt2, K());
}


template <class K>
inline 
typename K::FT
squared_distance_to_plane(const Vector_3<K> & normal,
			  const Vector_3<K> & diff)
{
  return internal::squared_distance_to_plane(normal, diff, K());
}


template <class K>
inline
typename K::FT
squared_distance_to_line(const Vector_3<K> & dir,
			 const Vector_3<K> & diff)
{
  return internal::squared_distance_to_line(dir, diff, K());
}



} //namespace CGAL


#endif
