// Copyright (c) 1998-2004  
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


#ifndef CGAL_DISTANCE_3_2_H
#define CGAL_DISTANCE_3_2_H

#include <CGAL/squared_distance_3_0.h>

#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Plane_3.h>

namespace CGAL {

namespace internal {

template <class K>
bool
contains_vector(const typename K::Plane_3 &pl, 
		const typename K::Vector_3 &vec,
		const K&)
{
  typedef typename K::RT RT;
  return pl.a()*vec.hx() + pl.b()*vec.hy() + pl.c() * vec.hz() == RT(0);
}


template <class K>
inline typename K::FT
squared_distance(
    const typename K::Point_3 & pt,
    const typename K::Plane_3 & plane,
    const K& k)
{
  typename K::Construct_vector_3 construct_vector;
  typedef typename K::Vector_3 Vector_3;
  Vector_3 diff = construct_vector(plane.point(), pt);
  return squared_distance_to_plane(plane.orthogonal_vector(), diff, k);
}



template <class K>
inline typename K::FT
squared_distance(
    const typename K::Plane_3 & plane,
    const typename K::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, plane, k);
}

template <class K>
typename K::FT
squared_distance(
    const typename K::Line_3 &line,
    const typename K::Plane_3 &plane,
    const K& k)
{
    typedef typename K::FT FT;
    if (contains_vector(plane, line.direction().vector(), k))
        return squared_distance(plane, line.point(), k);
    return FT(0);
}


template <class K>
inline typename K::FT
squared_distance(
    const typename K::Plane_3 & p,
    const typename K::Line_3 & line,
    const K& k)
{
    return squared_distance(line, p, k);
}

template <class K>
typename K::FT
squared_distance(
    const typename K::Ray_3 &ray,
    const typename K::Plane_3 &plane,
    const K& k)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 &start = ray.start();
    const Point_3 &planepoint = plane.point();
    Vector_3 start_min_pp = construct_vector(planepoint, start);
    Vector_3 end_min_pp = ray.direction().vector();
    const Vector_3 &normal = plane.orthogonal_vector();
    RT sdm_rs2pp = wdot(normal, start_min_pp, k);
    RT sdm_re2pp = wdot(normal, end_min_pp, k);
    switch (CGAL_NTS sign(sdm_rs2pp)) {
    case -1:
        if (sdm_re2pp > RT(0))
            return FT(0);
        return squared_distance_to_plane(normal, start_min_pp, k);
    case 0:
    default:
        return FT(0);
    case 1:
        if (sdm_re2pp < RT(0))
            return FT(0);
        return squared_distance_to_plane(normal, start_min_pp, k);
    }
}


template <class K>
inline typename K::FT
squared_distance(
    const typename K::Plane_3 & plane,
    const typename K::Ray_3 & ray,
    const K& k)
{
    return squared_distance(ray, plane, k);
}

template <class K>
typename K::FT
squared_distance(
    const typename K::Segment_3 &seg,
    const typename K::Plane_3 &plane,
    const K& k)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 &start = seg.start();
    const Point_3 &end = seg.end();
    if (start == end)
        return squared_distance(start, plane, k);
    const Point_3 &planepoint = plane.point();
    Vector_3 start_min_pp = construct_vector(planepoint, start);
    Vector_3 end_min_pp = construct_vector(planepoint, end);
    const Vector_3 &normal = plane.orthogonal_vector();
    RT sdm_ss2pp = wdot(normal, start_min_pp, k);
    RT sdm_se2pp = wdot(normal, end_min_pp, k);
    switch (CGAL_NTS sign(sdm_ss2pp)) {
    case -1:
        if (sdm_se2pp >= RT(0))
            return FT(0);
        if (sdm_ss2pp * end_min_pp.hw() >= sdm_se2pp * start_min_pp.hw())
            return squared_distance_to_plane(normal, start_min_pp, k);
        else
            return squared_distance_to_plane(normal, end_min_pp, k);
    case 0:
    default:
        return FT(0);
    case 1:
        if (sdm_se2pp <= RT(0))
            return FT(0);
        if (sdm_ss2pp  * end_min_pp.hw() <= sdm_se2pp * start_min_pp.hw())
            return squared_distance_to_plane(normal, start_min_pp, k);
        else
            return squared_distance_to_plane(normal, end_min_pp, k);
    }
}


template <class K>
inline typename K::FT
squared_distance(
    const typename K::Plane_3 & plane,
    const typename K::Segment_3 & seg,
    const K& k)
{
    return squared_distance(seg, plane, k);
}


} // namespace internal


template <class K>
bool
contains_vector(const Plane_3<K> &pl, const Vector_3<K> &vec)
{
  return internal::contains_vector(pl,vec, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Point_3<K> & pt,
    const Plane_3<K> & plane)
{
  return internal::squared_distance(pt, plane, K());
}



template <class K>
inline 
typename K::FT
squared_distance(
    const Plane_3<K> & plane,
    const Point_3<K> & pt)
{
    return internal::squared_distance(pt, plane, K());
}

template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> &line,
    const Plane_3<K> &plane)
{
    return internal::squared_distance(line, plane, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Plane_3<K> & p,
    const Line_3<K> & line)
{
    return internal::squared_distance(line, p, K());
}

template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> &ray,
    const Plane_3<K> &plane)
{
  return internal::squared_distance(ray, plane, K());
}



template <class K>
inline 
typename K::FT
squared_distance(
    const Plane_3<K> & plane,
    const Ray_3<K> & ray)
{
    return internal::squared_distance(ray, plane, K());
}

template <class K>
inline
typename K::FT
squared_distance(
    const Segment_3<K> &seg,
    const Plane_3<K> &plane)
{
  return internal::squared_distance(seg, plane, K());

}


template <class K>
inline 
typename K::FT
squared_distance(
    const Plane_3<K> & plane,
    const Segment_3<K> & seg)
{
    return internal::squared_distance(seg, plane, K());
}


} //namespace CGAL


#endif
