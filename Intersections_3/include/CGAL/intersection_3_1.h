// Copyright (c) 2010 
// GeometryFactory (France).
// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Sebastien Loriot <Sebastien.Loriot@geometryfactory.com>

#ifndef CGAL_INTERSECTION_3_1_H
#define CGAL_INTERSECTION_3_1_H

#include <CGAL/Object.h>
#include <CGAL/bbox_intersection_3.h>

namespace CGAL {

template <class R>
Object
intersection(const Plane_3<R> &plane1, const Plane_3<R>&plane2);

template <class R>
inline bool
do_intersect(const Plane_3<R> &plane1, const Plane_3<R>&plane2)
{
    return ! intersection(plane1, plane2).is_empty();
}


template <class R>
Object
intersection(const Plane_3<R> &plane1, const Plane_3<R>&plane2,
             const Plane_3<R>&plane3);

template <class R>
inline bool
do_intersect(const Plane_3<R> &plane1, const Plane_3<R>&plane2,
             const Plane_3<R>&plane3)
{
    return ! intersection(plane1, plane2, plane3).is_empty();
}


template <class R>
Object
intersection(const Plane_3<R> &plane, const Line_3<R>&line);

template <class R>
inline
Object
intersection(const Line_3<R>&line, const Plane_3<R> &plane)
{
    return intersection(plane,line);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p2, const Line_3<R> &p1);


template <class R>
inline bool
do_intersect(const Line_3<R> &p1, const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}



template <class R>
Object
intersection(const Plane_3<R> &plane, const Ray_3<R>&ray);

template <class R>
inline
Object
intersection(const Ray_3<R>&ray, const Plane_3<R> &plane)
{
    return intersection(plane,ray);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p1, const Ray_3<R> &p2);


template <class R>
inline bool
do_intersect(
    const Ray_3<R> &p1,
    const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}


template <class R>
Object
intersection(const Plane_3<R> &plane, const Segment_3<R>&seg);

template <class R>
inline
Object
intersection(const Segment_3<R>&seg, const Plane_3<R> &plane)
{
    return intersection(plane,seg);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p1, const Segment_3<R> &p2);


template <class R>
inline bool
do_intersect(const Segment_3<R> &p1,
	     const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}

template <class R>
Object
intersection(const Plane_3<R> &plane, const Triangle_3<R>&tri);

template <class R>
inline
Object
intersection(const Triangle_3<R>&tri, const Plane_3<R> &plane)
{
    return intersection(plane,tri);
}

template <class R>
Object
intersection(const Line_3<R> &line,
	     const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
	     const Line_3<R> &line)
{
    return intersection(line, box);
}



template <class R>
Object
intersection(const Ray_3<R> &ray,
	     const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
	     const Ray_3<R> &ray)
{
    return intersection(ray, box);
}



template <class R>
Object
intersection(const Segment_3<R> &seg,
	     const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
	     const Segment_3<R> &seg)
{
    return intersection(seg, box);
}


template <class R>
Object
intersection(const Line_3<R> &line,
	     const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
	     const Line_3<R> &line)
{
    return intersection(line, box);
}



template <class R>
Object
intersection(const Ray_3<R> &ray,
	     const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
	     const Ray_3<R> &ray)
{
    return intersection(ray, box);
}



template <class R>
Object
intersection(const Segment_3<R> &seg,
	     const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
	     const Segment_3<R> &seg)
{
    return intersection(seg, box);
}



template <class R>
Object
intersection(const Iso_cuboid_3<R> &box1,
	     const Iso_cuboid_3<R> &box2) ;

template <class R>
Object
intersection(const Line_3<R> &l1,
             const Line_3<R> &l2);

template <class R>
Object
intersection(const Segment_3<R> &,
             const Segment_3<R> &);

template <class R>
Object
intersection(const Line_3<R> &,
             const Segment_3<R> &);
             
template <class R>
Object
intersection(const Segment_3<R> &s,
             const Line_3<R> &l)
{return intersection(l,s);}

template <class R>
Object
intersection(const Line_3<R> &,
             const Ray_3<R> &);
             
template <class R>
Object
intersection(const Ray_3<R> &r,
             const Line_3<R> &l)
{return intersection(l,r);}

template <class R>
Object
intersection(const Ray_3<R> &,
             const Segment_3<R> &);
             
template <class R>
Object
intersection(const Segment_3<R> &s,
             const Ray_3<R> &r)
{return intersection(r,s);}

template <class R>
Object
intersection(const Ray_3<R> &,
             const Ray_3<R> &);

template <class R>
Object
intersection(const Sphere_3<R> &s1,
             const Sphere_3<R> &s2);


template <class R>
Object
intersection(const Plane_3<R> &p,
             const Sphere_3<R> &s);

template <class R>
inline Object
intersection(const Sphere_3<R> &s,
             const Plane_3<R> &p) {
  return intersection(p, s);
}

template <class R>
bool
do_intersect(const Line_3<R> &l1,
             const Line_3<R> &l2);

template <class R>
bool
do_intersect(const Segment_3<R> &,
             const Segment_3<R> &);

template <class R>
bool
do_intersect(const Line_3<R> &,
             const Segment_3<R> &);
             
template <class R>
bool
do_intersect(const Segment_3<R> &s,
             const Line_3<R> &l)
{return do_intersect(l,s);}

template <class R>
bool
do_intersect(const Line_3<R> &,
             const Ray_3<R> &);
             
template <class R>
bool
do_intersect(const Ray_3<R> &r,
             const Line_3<R> &l)
{return do_intersect(l,r);}

template <class R>
bool
do_intersect(const Ray_3<R> &,
             const Segment_3<R> &);
             
template <class R>
bool
do_intersect(const Segment_3<R> &s,
             const Ray_3<R> &r)
{return do_intersect(r,s);}

template <class R>
bool
do_intersect(const Ray_3<R> &r1,
             const Ray_3<R> &r2);

template <class R>
bool
do_intersect(const Sphere_3<R> &s1,
             const Sphere_3<R> &s2);

template <class R>
bool
do_intersect(const Plane_3<R> &p,
             const Sphere_3<R> &s);

template <class R>
inline bool
do_intersect(const Sphere_3<R> &s,
             const Plane_3<R> &p) {
  return do_intersect(p, s);
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &i, const Iso_cuboid_3<R> &j)
{
	return ! CGAL::intersection(i, j).empty();
}

template <class R>
inline bool
do_intersect(const Line_3<R> &l, const Iso_cuboid_3<R> &j)
{
	return ! CGAL::intersection(l, j).empty();
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Line_3<R> &l)
{
	return ! CGAL::intersection(l, j).empty();
}

template <class R>
inline bool
do_intersect(const Ray_3<R> &r, const Iso_cuboid_3<R> &j)
{
	return ! CGAL::intersection(r, j).empty();
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Ray_3<R> &r)
{
	return ! CGAL::intersection(r, j).empty();
}

template <class R>
inline bool
do_intersect(const Segment_3<R> &s, const Iso_cuboid_3<R> &j)
{
	return ! CGAL::intersection(s, j).empty();
}

template <class R>
inline bool
do_intersect(const Iso_cuboid_3<R> &j, const Segment_3<R> &s)
{
	return ! CGAL::intersection(s, j).empty();
}

} //namespace CGAL

#include <CGAL/Intersections_3/intersection_3_1_impl.h>

#endif // CGAL_INTERSECTION_3_1_H
