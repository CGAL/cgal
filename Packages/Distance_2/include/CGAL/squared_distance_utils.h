// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/squared_distance_utils.h
// source        : sqdistance_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_SQUARED_DISTANCE_UTILS_H
#define CGAL_SQUARED_DISTANCE_UTILS_H

#include <CGAL/determinant.h>
#include <CGAL/wmult.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

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
typename K::RT wdot(const typename K::Point_2 &p,
		    const typename K::Point_2 &q,
		    const typename K::Point_2 &r,
		    const K&)
{
    K* pR = 0;
    return  (wmult(pR, p.hx(),q.hw()) - wmult(pR, q.hx(),p.hw()))
          * (wmult(pR, r.hx(),q.hw()) - wmult(pR, q.hx(),r.hw()))
          + (wmult(pR, p.hy(),q.hw()) - wmult(pR, q.hy(),p.hw()))
          * (wmult(pR, r.hy(),q.hw()) - wmult(pR, q.hy(),r.hw()));
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
    return det3x3_by_formula(
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
Orientation orientation(const typename K::Vector_2 &u,
			const typename K::Vector_2 &v,
			const K& k)
{
    typedef typename K::RT RT;
    RT wcr = wcross(u,v, k);
    return (wcr > RT(0)) ? COUNTERCLOCKWISE :
           (wcr < RT(0)) ? CLOCKWISE
                            : COLLINEAR;
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

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_SQUARED_DISTANCE_UTILS_H
