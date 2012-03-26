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


#ifndef CGAL_DISTANCE_3_1_H
#define CGAL_DISTANCE_3_1_H


#include <CGAL/squared_distance_3_0.h>

#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>


namespace CGAL {

namespace internal {

template <class K>
typename K::FT
squared_distance(
    const typename K::Point_3 &pt,
    const typename K::Line_3 &line,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typename K::Construct_vector_3 construct_vector;
  Vector_3 dir(line.direction().vector());
  Vector_3 diff = construct_vector(line.point(), pt);
  return internal::squared_distance_to_line(dir, diff, k);
}

template <class K>
inline
typename K::FT
squared_distance(
    const typename K::Line_3 & line,
    const typename K::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, line, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename K::Point_3 &pt,
    const typename K::Ray_3 &ray,
    const K& k)
{
  typename K::Construct_vector_3 construct_vector;
  typedef typename K::Vector_3 Vector_3;

    Vector_3 diff = construct_vector(ray.source(), pt);
    const Vector_3 &dir = ray.direction().vector();
    if (!is_acute_angle(dir,diff, k) )
        return (typename K::FT)(diff*diff);
    return squared_distance_to_line(dir, diff, k);
}


template <class K>
inline
typename K::FT
squared_distance(
    const typename K::Ray_3 & ray,
    const typename K::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, ray, k);
}




template <class K>
typename K::FT
squared_distance(
    const typename K::Point_3 &pt,
    const typename K::Segment_3 &seg,
    const K& k,
    const Homogeneous_tag)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    // assert that the segment is valid (non zero length).
    Vector_3 diff = construct_vector(seg.source(), pt);
    Vector_3 segvec = construct_vector(seg.source(), seg.target());
    RT d = wdot(diff,segvec, k);
    if (d <= (RT)0)
        return (FT(diff*diff));
    RT e = wdot(segvec,segvec, k);
    if ( (d * segvec.hw()) > (e * diff.hw()))
        return squared_distance(pt, seg.target(), k);

    Vector_3 wcr = wcross(segvec, diff, k);
    return FT(wcr*wcr)/FT(e * diff.hw() * diff.hw());
}

template <class K>
typename K::FT
squared_distance(
    const typename K::Point_3 &pt,
    const typename K::Segment_3 &seg,
    const K& k,
    const Cartesian_tag&)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    // assert that the segment is valid (non zero length).
    Vector_3 diff = construct_vector(seg.source(), pt);
    Vector_3 segvec = construct_vector(seg.source(), seg.target());
    RT d = wdot(diff,segvec, k);
    if (d <= (RT)0)
        return (FT(diff*diff));
    RT e = wdot(segvec,segvec, k);
    if (d > e)
        return squared_distance(pt, seg.target(), k);

    Vector_3 wcr = wcross(segvec, diff, k);
    return FT(wcr*wcr)/e;
}


template <class K>
inline
typename K::FT
squared_distance(
    const typename K::Point_3 &pt,
    const typename K::Segment_3 &seg,
    const K& k)
{ 
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return squared_distance(pt, seg, k, tag);
}


template <class K>
inline 
typename K::FT
squared_distance(
    const typename K::Segment_3 & seg,
    const typename K::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, seg, k);
}




template <class K>
typename K::FT
squared_distance_parallel(
    const typename K::Segment_3 &seg1,
    const typename K::Segment_3 &seg2,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
    const Vector_3 &dir1 = seg1.direction().vector();
    const Vector_3 &dir2 = seg2.direction().vector();
 
    if (same_direction(dir1, dir2, k)) {
        if (!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
            return squared_distance(seg1.target(), seg2.source(), k);
        if (!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
            return squared_distance(seg1.source(), seg2.target(), k);
    } else {
        if (!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
            return squared_distance(seg1.target(), seg2.target(), k);
        if (!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
            return squared_distance(seg1.source(), seg2.source(), k);
    }
    return squared_distance(seg2.source(), seg1.supporting_line(), k);
}



template <class K>
inline
typename K::RT 
_distance_measure_sub(typename K::RT startwdist, typename K::RT endwdist,
			 const typename K::Vector_3 &start, 
			 const typename K::Vector_3 &end,
			 const K&)
{
    return  CGAL_NTS abs(wmult((K*)0, startwdist, end.hw())) -
            CGAL_NTS abs(wmult((K*)0, endwdist, start.hw()));
}


template <class K>
typename K::FT
squared_distance(
    const typename K::Segment_3 &seg1,
    const typename K::Segment_3 &seg2,
    const K& k)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 &start1 = seg1.source();
    const Point_3 &start2 = seg2.source();
    const Point_3 &end1 = seg1.target();
    const Point_3 &end2 = seg2.target();

    if (start1 == end1)
        return squared_distance(start1, seg2, k);
    if (start2 == end2)
        return squared_distance(start2, seg1, k);
    
    Vector_3 dir1, dir2, normal;
    dir1 = seg1.direction().vector();
    dir2 = seg2.direction().vector();
    normal = wcross(dir1, dir2, k);
    if (is_null(normal, k))
        return squared_distance_parallel(seg1, seg2, k);
    
    bool crossing1, crossing2;
    RT sdm_s1to2, sdm_e1to2, sdm_s2to1, sdm_e2to1;
    Vector_3 perpend1, perpend2, s2mins1, e2mins1, e1mins2;
    perpend1 = wcross(dir1, normal, k);
    perpend2 = wcross(dir2, normal, k);
    s2mins1 = construct_vector(start1, start2);
    e2mins1 = construct_vector(start1, end2);
    e1mins2 = construct_vector(start2, end1);
    sdm_s1to2 = -RT(wdot(perpend2, s2mins1, k));
    sdm_e1to2 = wdot(perpend2, e1mins2, k);
    sdm_s2to1 = wdot(perpend1, s2mins1, k);
    sdm_e2to1 = wdot(perpend1, e2mins1, k);
    
    if (sdm_s1to2 < RT(0)) {
        crossing1 = (sdm_e1to2 >= RT(0));
    } else {
        if (sdm_e1to2 <= RT(0)) {
            crossing1 = true;
        } else {
            crossing1 = (sdm_s1to2 == RT(0));
        }
    }
    if (sdm_s2to1 < RT(0)) {
        crossing2 = (sdm_e2to1 >= RT(0));
    } else {
        if (sdm_e2to1 <= RT(0)) {
            crossing2 = true;
        } else {
            crossing2 = (sdm_s2to1 == RT(0));
        }
    }
    
    if (crossing1) {
        if (crossing2) {
            return squared_distance_to_plane(normal, s2mins1, k);
        }
    
        RT dm;
        dm = _distance_measure_sub(
                  sdm_s2to1, sdm_e2to1, s2mins1, e2mins1, k);
        if (dm < RT(0)) {
            return squared_distance(start2, seg1, k);
        } else {
            if (dm > RT(0)) {
                return squared_distance(end2, seg1, k);
            } else {
                // should not happen with exact arithmetic.
                return squared_distance_parallel(seg1, seg2, k);
            }
        }
    } else {
        if (crossing2) {
            RT dm;
            dm =_distance_measure_sub(
                 sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
            if (dm < RT(0)) {
                return squared_distance(start1, seg2, k);
            } else {
                if (dm > RT(0)) {
                    return squared_distance(end1, seg2, k);
                } else {
                    // should not happen with exact arithmetic.
                    return squared_distance_parallel(seg1, seg2, k);
                }
            }
        } else {
            FT min1, min2;
            RT dm;
            dm = _distance_measure_sub(
                     sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
            if (dm == RT(0)) // should not happen with exact arithmetic.
               return squared_distance_parallel(seg1, seg2, k);
            min1 = (dm < RT(0)) ?
                squared_distance(seg1.source(), seg2, k):
                squared_distance(end1, seg2, k);
            dm = _distance_measure_sub(
                     sdm_s2to1, sdm_e2to1, s2mins1, e2mins1, k);
            if (dm == RT(0)) // should not happen with exact arithmetic.
                return squared_distance_parallel(seg1, seg2, k);
            min2 = (dm < RT(0)) ?
                squared_distance(start2, seg1, k):
                squared_distance(end2, seg1, k);
            return (min1 < min2) ? min1 : min2;
        }
    }
    
}






template <class K>
typename K::FT
squared_distance_parallel(
    const typename K::Segment_3 &seg,
    const typename K::Ray_3 &ray,
    const K& k)
{

  typedef typename K::Vector_3 Vector_3;
  bool same_direction;
  const Vector_3 &dir1 = seg.direction().vector();
  const Vector_3 &dir2 = ray.direction().vector();
  if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
    same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
  } else {
    same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
  }
  if (same_direction) {
    if (!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
      return squared_distance(seg.target(), ray.source(), k);
  } else {
    if (!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
      return squared_distance(seg.source(), ray.source(), k);
  }
  return squared_distance(ray.source(), seg.supporting_line(), k);
}


template <class K>
typename K::FT
squared_distance(
    const typename K::Segment_3 &seg,
    const typename K::Ray_3 &ray,
    const K& k)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Point_3 Point_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 & ss = seg.source();
    const Point_3 & se = seg.target();
    if (ss == se)
        return squared_distance(ss, ray, k);
    Vector_3 raydir, segdir, normal;
    raydir = ray.direction().vector();
    segdir = seg.direction().vector();
    normal = wcross(segdir, raydir, k);
    if (is_null(normal, k))
        return squared_distance_parallel(seg, ray, k);

    bool crossing1, crossing2;
    RT sdm_ss2r, sdm_se2r, sdm_rs2s, sdm_re2s;
    Vector_3 perpend2seg, perpend2ray, ss_min_rs, se_min_rs;
    perpend2seg = wcross(segdir, normal, k);
    perpend2ray = wcross(raydir, normal, k);
    ss_min_rs = construct_vector(ray.source(), ss);
    se_min_rs = construct_vector(ray.source(), se);
    sdm_ss2r = wdot(perpend2ray, ss_min_rs, k);
    sdm_se2r = wdot(perpend2ray, se_min_rs, k);
    if (sdm_ss2r < RT(0)) {
        crossing1 = (sdm_se2r >= RT(0));
    } else {
        if (sdm_se2r <= RT(0)) {
            crossing1 = true;
        } else {
            crossing1 = (sdm_ss2r == RT(0));
        }
    }

    sdm_rs2s = -RT(wdot(perpend2seg, ss_min_rs, k));
    sdm_re2s = wdot(perpend2seg, raydir, k);
    if (sdm_rs2s < RT(0)) {
        crossing2 = (sdm_re2s >= RT(0));
    } else {
        if (sdm_re2s <= RT(0)) {
            crossing2 = true;
        } else {
            crossing2 = (sdm_rs2s == RT(0));
        }
    }

    if (crossing1) {
        if (crossing2) {
            return squared_distance_to_plane(normal, ss_min_rs, k);
        }
        return squared_distance(ray.source(), seg, k);
    } else {
        if (crossing2) {
            RT dm;
            dm = _distance_measure_sub(
                    sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
            if (dm < RT(0)) {
                return squared_distance(ss, ray, k);
            } else {
                if (dm > RT(0)) {
                    return squared_distance(se, ray, k);
                } else {
                    // parallel, should not happen (no crossing)
                    return squared_distance_parallel(seg, ray, k);
                }
            }
        } else {
            FT min1, min2;
            RT dm;
            dm = _distance_measure_sub(
                    sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
            if (dm == RT(0))
                return squared_distance_parallel(seg, ray, k);
            min1 = (dm < RT(0))
                 ? squared_distance(ss, ray, k)
                 : squared_distance(se, ray, k);
            min2 = squared_distance(ray.source(), seg, k);
            return (min1 < min2) ? min1 : min2;
        }
    }
}



template <class K>
inline
typename K::FT
squared_distance(
    const typename K::Ray_3 & ray,
    const typename K::Segment_3 & seg,
    const K& k)
{
    return squared_distance(seg, ray, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename K::Segment_3 &seg,
    const typename K::Line_3 &line,
    const K& k)
{
    typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    const Point_3 &linepoint = line.point();
    const Point_3 &start = seg.source();
    const Point_3 &end = seg.target();

    if (start == end)
        return squared_distance(start, line, k);
    Vector_3 linedir = line.direction().vector();
    Vector_3 segdir = seg.direction().vector();
    Vector_3 normal = wcross(segdir, linedir, k);
    if (is_null(normal, k))
        return squared_distance_to_line(linedir, 
					construct_vector(linepoint,start), k);

    bool crossing;
    RT sdm_ss2l, sdm_se2l;
    Vector_3 perpend2line, start_min_lp, end_min_lp;
    perpend2line = wcross(linedir, normal, k);
    start_min_lp = construct_vector(linepoint, start);
    end_min_lp = construct_vector(linepoint, end);
    sdm_ss2l = wdot(perpend2line, start_min_lp, k);
    sdm_se2l = wdot(perpend2line, end_min_lp, k);
    if (sdm_ss2l < RT(0)) {
        crossing = (sdm_se2l >= RT(0));
    } else {
        if (sdm_se2l <= RT(0)) {
            crossing = true;
        } else {
            crossing = (sdm_ss2l == RT(0));
        }
    }

    if (crossing) {
        return squared_distance_to_plane(normal, start_min_lp, k);
    } else {
        RT dm;
        dm = _distance_measure_sub(
                sdm_ss2l, sdm_se2l, start_min_lp, end_min_lp, k);
        if (dm <= RT(0)) {
            return squared_distance_to_line(linedir, start_min_lp, k);
        } else {
            return squared_distance_to_line(linedir, end_min_lp, k);
        }
    }
}


template <class K>
inline 
typename K::FT
squared_distance(
    const typename K::Line_3 & line,
    const typename K::Segment_3 & seg,
    const K& k)
{
    return squared_distance(seg, line, k);
}




template <class K>
typename K::FT
ray_ray_squared_distance_parallel(
    const typename K::Vector_3 &ray1dir,
    const typename K::Vector_3 &ray2dir,
    const typename K::Vector_3 &s1_min_s2,
    const K& k)
{
  if (!is_acute_angle(ray2dir, s1_min_s2, k)) {
    if (!same_direction(ray1dir, ray2dir, k))
      return (typename K::FT)(s1_min_s2*s1_min_s2);
  }
  return squared_distance_to_line(ray1dir, s1_min_s2, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename K::Ray_3 &ray1,
    const typename K::Ray_3 &ray2,
    const K& k)
{
  typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 & s1 = ray1.source();
    const Point_3 & s2 = ray2.source();
    Vector_3 dir1, dir2, normal;
    dir1 = ray1.direction().vector();
    dir2 = ray2.direction().vector();
    normal = wcross(dir1, dir2, k);
    Vector_3 s1_min_s2 = construct_vector(s2, s1);
    if (is_null(normal, k))
        return ray_ray_squared_distance_parallel(dir1, dir2, s1_min_s2, k);

    bool crossing1, crossing2;
    RT sdm_s1_2, sdm_s2_1;
    Vector_3 perpend1, perpend2;
    perpend1 = wcross(dir1, normal, k);
    perpend2 = wcross(dir2, normal, k);

    sdm_s1_2 = wdot(perpend2, s1_min_s2, k);
    if (sdm_s1_2 < RT(0)) {
        crossing1 = (RT(wdot(perpend2, dir1, k)) >= RT(0));
    } else {
        if (RT(wdot(perpend2, dir1, k)) <= RT(0)) {
            crossing1 = true;
        } else {
            crossing1 = (sdm_s1_2 == RT(0));
        }
    }
    sdm_s2_1 = -RT(wdot(perpend1, s1_min_s2, k));
    if (sdm_s2_1 < RT(0)) {
        crossing2 = (RT(wdot(perpend1, dir2, k)) >= RT(0));
    } else {
        if (RT(wdot(perpend1, dir2, k)) <= RT(0)) {
            crossing2 = true;
        } else {
            crossing2 = (sdm_s2_1 == RT(0));
        }
    }
    if (crossing1) {
        if (crossing2)
            return squared_distance_to_plane(normal, s1_min_s2, k);
        return squared_distance(ray2.source(), ray1, k);
    } else {
        if (crossing2) {
            return squared_distance(ray1.source(), ray2, k);
        } else {
          FT min1, min2;
            min1 = squared_distance(ray1.source(), ray2, k);
            min2 = squared_distance(ray2.source(), ray1, k);
            return (min1 < min2) ? min1 : min2;
        }
    }
}





template <class K>
typename K::FT
squared_distance(
    const typename K::Line_3 &line,
    const typename K::Ray_3 &ray,
    const K& k)
{
  typename K::Construct_vector_3 construct_vector;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    const Point_3 & rs =ray.source();
    Vector_3 raydir, linedir, normal;
    linedir = line.direction().vector();
    raydir = ray.direction().vector();
    normal = wcross(raydir, linedir, k);
    Vector_3 rs_min_lp = construct_vector(line.point(), rs);
    if (is_null(normal, k))
        return squared_distance_to_line(linedir, rs_min_lp, k);

    bool crossing;
    RT sdm_sr_l;
    Vector_3 perpend2l;
    perpend2l = wcross(linedir, normal, k);

    sdm_sr_l = wdot(perpend2l, rs_min_lp, k);
    if (sdm_sr_l < RT(0)) {
        crossing = (RT(wdot(perpend2l, raydir, k)) >= RT(0));
    } else {
        if (RT(wdot(perpend2l, raydir, k)) <= RT(0)) {
            crossing = true;
        } else {
            crossing = (sdm_sr_l == RT(0));
        }
    }

    if (crossing)
        return squared_distance_to_plane(normal, rs_min_lp, k);
    else
        return squared_distance_to_line(linedir, rs_min_lp, k);
}


template <class K>
inline typename K::FT
squared_distance(
    const typename K::Ray_3 & ray,
    const typename K::Line_3 & line,
    const K& k)
{
    return squared_distance(line, ray, k);
}




template <class K>
typename K::FT
squared_distance(
    const typename K::Line_3 &line1,
    const typename K::Line_3 &line2,
    const K& k)
{   
  typename K::Construct_vector_3 construct_vector;
    typedef typename K::Vector_3 Vector_3;

    Vector_3 dir1, dir2, normal, diff;
    dir1 = line1.direction().vector();
    dir2 = line2.direction().vector();
    normal = wcross(dir1, dir2, k);
    diff = construct_vector(line1.point(), line2.point());
    if (is_null(normal, k))
        return squared_distance_to_line(dir2, diff, k);
    return squared_distance_to_plane(normal, diff, k);
}



} // namespace internal



template <class K>
inline
typename K::FT
squared_distance(const Point_3<K> &pt,
		 const Line_3<K> &line)
{
  return internal::squared_distance(pt, line, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> & line,
    const Point_3<K> & pt)
{
  return internal::squared_distance(pt, line, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Point_3<K> &pt,
    const Ray_3<K> &ray)
{
    return internal::squared_distance(pt, ray, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Point_3<K> & pt)
{
    return internal::squared_distance(pt, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Point_3<K> &pt,
    const Segment_3<K> &seg)
{
  return internal::squared_distance(pt, seg, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Segment_3<K> & seg,
    const Point_3<K> & pt)
{
    return internal::squared_distance(pt, seg, K());
}




template <class K>
inline
typename K::FT
squared_distance_parallel(
    const Segment_3<K> &seg1,
    const Segment_3<K> &seg2)
{
  return internal::squared_distance_parallel(seg1, seg2, K());
}




template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K> &seg1,
		 const Segment_3<K> &seg2)
{
  return internal::squared_distance(seg1, seg2, K());
}






template <class K>
inline
typename K::FT
squared_distance_parallel(
    const Segment_3<K> &seg,
    const Ray_3<K> &ray)
{
  return internal::squared_distance_parallel(ray,seg, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Segment_3<K> &seg,
    const Ray_3<K> &ray)
{
  return internal::squared_distance(ray, seg, K());
}



template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Segment_3<K> & seg)
{
    return internal::squared_distance(seg, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Segment_3<K> &seg,
    const Line_3<K> &line)
{
  return internal::squared_distance(seg, line, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Line_3<K> & line,
    const Segment_3<K> & seg)
{
    return internal::squared_distance(seg, line, K());
}




template <class K>
inline
typename K::FT
ray_ray_squared_distance_parallel(
    const Vector_3<K> &ray1dir,
    const Vector_3<K> &ray2dir,
    const Vector_3<K> &s1_min_s2)
{
  return internal::ray_ray_squared_distance_parallel(ray1dir, ray2dir, 
						  s1_min_s2, K());
}

template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> &ray1,
    const Ray_3<K> &ray2)
{
  return internal::squared_distance(ray1, ray2, K());
}





template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> &line,
    const Ray_3<K> &ray)
{
  return internal::squared_distance(line, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Line_3<K> & line)
{
    return internal::squared_distance(line, ray, K());
}




template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> &line1,
    const Line_3<K> &line2)
{  
    return internal::squared_distance(line1, line2, K());
}



} //namespace CGAL


#endif
