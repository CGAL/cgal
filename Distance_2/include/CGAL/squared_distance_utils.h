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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_SQUARED_DISTANCE_UTILS_H
#define CGAL_SQUARED_DISTANCE_UTILS_H

#include <CGAL/determinant.h>
#include <CGAL/wmult.h>

namespace CGAL {

namespace internal {

template <class K>
bool is_null(const  typename K::Vector_2 &v, const K&)
{
    typedef typename K::RT RT;
    return v.hx()==RT(0) && v.hy()==RT(0);
}


template <class K>
typename K::RT
wdot(const typename K::Vector_2 &u, 
     const typename K::Vector_2 &v,
     const K&)
{
    return  (u.hx()*v.hx() + u.hy()*v.hy());
}



template <class K>
typename K::RT wdot_tag(const typename K::Point_2 &p,
			const typename K::Point_2 &q,
			const typename K::Point_2 &r,
			const K&,
			const Cartesian_tag&)
{
  return  (p.x() - q.x()) * (r.x() - q.x())
          + (p.y() - q.y()) * (r.y() - q.y());
}


template <class K>
typename K::RT wdot_tag(const typename K::Point_2 &p,
			const typename K::Point_2 &q,
			const typename K::Point_2 &r,
			const K&,
			const Homogeneous_tag&)
{
  return  (p.hx() * q.hw() - q.hx() * p.hw())
          * (r.hx() * q.hw() - q.hx() * r.hw())
          + (p.hy() * q.hw() - q.hy() * p.hw())
            * (r.hy() * q.hw() - q.hy() * r.hw());
}


template <class K>
typename K::RT wdot(const typename K::Point_2 &p,
		    const typename K::Point_2 &q,
		    const typename K::Point_2 &r,
		    const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return wdot_tag(p, q, r, k, tag);
}



template <class K>
typename K::RT
wcross(const typename K::Vector_2 &u,
       const typename K::Vector_2 &v,
       const K&)
{
    return (typename K::RT)(u.hx()*v.hy() - u.hy()*v.hx());
}



template <class K>
inline
typename K::RT 
wcross_tag(const typename K::Point_2 &p,
	   const typename K::Point_2 &q,
	   const typename K::Point_2 &r,
	   const K&,
	   const Homogeneous_tag&)
{
    return CGAL::determinant(
        p.hx(), q.hx(), r.hx(),
        p.hy(), q.hy(), r.hy(),
        p.hw(), q.hw(), r.hw());
}



template <class K>
inline
typename K::FT 
wcross_tag(const typename K::Point_2 &p,
	   const typename K::Point_2 &q,
	   const typename K::Point_2 &r,
	   const K&,
	   const Cartesian_tag&)
{
  return (q.x()-p.x())*(r.y()-q.y()) - (q.y()-p.y())*(r.x()-q.x());
}


template <class K>
typename K::RT wcross(const typename K::Point_2 &p,
		      const typename K::Point_2 &q,
		      const typename K::Point_2 &r,
		      const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return wcross_tag(p, q, r, k, tag);

}



template <class K>
inline bool is_acute_angle(const typename K::Vector_2 &u,
			   const typename K::Vector_2 &v,
			   const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) > RT(0) ;
}

template <class K>
inline bool is_straight_angle(const typename K::Vector_2 &u,
			      const typename K::Vector_2 &v,
			      const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) == RT(0) ;
}

template <class K>
inline bool is_obtuse_angle(const typename K::Vector_2 &u,
			    const typename K::Vector_2 &v,
			    const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(u, v, k)) < RT(0) ;
}

template <class K>
inline bool is_acute_angle(const typename K::Point_2 &p,
			   const typename K::Point_2 &q, 
			   const typename K::Point_2 &r,
			   const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) > RT(0) ;
}

template <class K>
inline bool is_straight_angle(const typename K::Point_2 &p,
			      const typename K::Point_2 &q, 
			      const typename K::Point_2 &r,
			      const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) == RT(0) ;
}

template <class K>
inline bool is_obtuse_angle(const typename K::Point_2 &p,
			    const typename K::Point_2 &q, 
			    const typename K::Point_2 &r,
			    const K& k)
{
    typedef typename K::RT RT;
    return RT(wdot(p, q, r, k)) < RT(0) ;
}

template <class K>
inline bool counterclockwise(const typename K::Vector_2 &u,
			     const typename K::Vector_2 &v,
			     const K& k)
{
    typedef typename K::RT RT;
    return RT(wcross(u,v, k)) > RT(0);
}

template <class K>
inline bool left_turn(const typename K::Vector_2 &u,
		      const typename K::Vector_2 &v,
		      const K& k)
{
    typedef typename K::RT RT;
    return RT(wcross(u,v, k)) > RT(0);
}

template <class K>
inline bool clockwise(const typename K::Vector_2 &u,
		      const typename K::Vector_2 &v,
		      const K& k)
{
    typedef typename K::RT RT;
    return RT(wcross(u,v, k)) < RT(0);
}

template <class K>
inline bool right_turn(const typename K::Vector_2 &u,
		       const typename K::Vector_2 &v,
		       const K& k)
{
    typedef typename K::RT RT;
    return RT(wcross(u,v, k)) < RT(0);
}

template <class K>
inline bool collinear(const typename K::Vector_2 &u,
		      const typename K::Vector_2 &v,
		      const K& k)
{
    typedef typename K::RT RT;
    return RT(wcross(u,v, k)) == RT(0);
}

/*
the ordertype, right_turn, left_turn and collinear routines for points are
defined elsewhere.
*/
template <class K>
inline
bool
same_direction_tag(const typename K::Vector_2 &u,
		   const typename K::Vector_2 &v,
		   const K&,
		   const Cartesian_tag&)
{ 
  typedef typename K::FT FT;
  const FT& ux = u.x();
  const FT& uy = u.y();
   if (CGAL_NTS abs(ux) > CGAL_NTS abs(uy)) {
      return CGAL_NTS sign(ux) == CGAL_NTS sign(v.x());
  } else {
    return CGAL_NTS sign(uy) == CGAL_NTS sign(v.y());
  } 
}


template <class K>
inline
bool
same_direction_tag(const typename K::Vector_2 &u,
		   const typename K::Vector_2 &v,
		   const K&,
		   const Homogeneous_tag&)
{   
  typedef typename K::RT RT;
  const RT& uhx = u.hx();
  const RT& uhy = u.hy();
  if (CGAL_NTS abs(uhx) > CGAL_NTS abs(uhy)) {
      return CGAL_NTS sign(uhx) == CGAL_NTS sign(v.hx());
  } else {
    return CGAL_NTS sign(uhy) == CGAL_NTS sign(v.hy());
  }
}


template <class K>
inline
bool
same_direction(const typename K::Vector_2 &u,
	       const typename K::Vector_2 &v,
	       const K& k)
{  
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return same_direction_tag(u,v, k, tag);
}


} // namespace internal

} //namespace CGAL

#endif // CGAL_SQUARED_DISTANCE_UTILS_H
