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
// file          : include/CGAL/squared_distance_3_1.h
// source        : sqdistance_3.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_DISTANCE_3_1_H
#define CGAL_DISTANCE_3_1_H

#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>

#include <CGAL/utils.h>
#include <CGAL/Point_3.h>
#include <CGAL/enum.h>
#include <CGAL/wmult.h>
#include <CGAL/squared_distance_3_0.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Point_3 &pt,
    const typename CGAL_WRAP(K)::Line_3 &line,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  Vector_3 dir(line.direction().vector());
  Vector_3 diff = pt - line.point();
  return CGALi::squared_distance_to_line(dir, diff, k);
}

template <class K>
inline
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Line_3 & line,
    const typename CGAL_WRAP(K)::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, line, k);
}


template <class K>
extern
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Point_3 &pt,
    const typename CGAL_WRAP(K)::Ray_3 &ray,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;

    Vector_3 diff = pt-ray.start();
    const Vector_3 &dir = ray.direction().vector();
    if (!is_acute_angle(dir,diff, k) )
        return (typename K::FT)(diff*diff);
    return squared_distance_to_line(dir, diff, k);
}


template <class K>
inline
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Ray_3 & ray,
    const typename CGAL_WRAP(K)::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, ray, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Point_3 &pt,
    const typename CGAL_WRAP(K)::Segment_3 &seg,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    // assert that the segment is valid (non zero length).
    Vector_3 diff = pt-seg.start();
    Vector_3 segvec = seg.end()-seg.start();
    RT d = wdot(diff,segvec, k);
    if (d <= (RT)0)
        return (FT(diff*diff));
    RT e = wdot(segvec,segvec, k);
    if (wmult((K*)0 ,d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
        return squared_distance(pt, seg.end(), k);
    return squared_distance_to_line(segvec, diff, k);
}


template <class K>
inline 
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Segment_3 & seg,
    const typename CGAL_WRAP(K)::Point_3 & pt,
    const K& k)
{
    return squared_distance(pt, seg, k);
}




template <class K>
typename K::FT
squared_distance_parallel(
    const typename CGAL_WRAP(K)::Segment_3 &seg1,
    const typename CGAL_WRAP(K)::Segment_3 &seg2,
    const K k)
{
  typedef typename K::Vector_3 Vector_3;
    bool same_direction;
    const Vector_3 &dir1 = seg1.direction().vector();
    const Vector_3 &dir2 = seg2.direction().vector();
    if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
        if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hz())) {
            same_direction =
                (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
        } else {
            same_direction =
                (CGAL_NTS sign(dir1.hz()) == CGAL_NTS sign(dir2.hz()));
        }
    } else {
        if (CGAL_NTS abs(dir1.hy()) > CGAL_NTS abs(dir1.hz())) {
            same_direction =
                (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
        } else {
            same_direction =
                (CGAL_NTS sign(dir1.hz()) == CGAL_NTS sign(dir2.hz()));
        }
    }
    if (same_direction) {
        if (!is_acute_angle(seg1.start(), seg1.end(), seg2.start(), k))
            return squared_distance(seg1.end(), seg2.start(), k);
        if (!is_acute_angle(seg1.end(), seg1.start(), seg2.end(), k))
            return squared_distance(seg1.start(), seg2.end(), k);
    } else {
        if (!is_acute_angle(seg1.start(), seg1.end(), seg2.end(), k))
            return squared_distance(seg1.end(), seg2.end(), k);
        if (!is_acute_angle(seg1.end(), seg1.start(), seg2.start(), k))
            return squared_distance(seg1.start(), seg2.start(), k);
    }
    return squared_distance(seg2.start(), seg1.supporting_line(), k);
}



template <class RT, class K>
inline
RT 
_distance_measure_sub(RT startwdist, RT endwdist,
			 const typename CGAL_WRAP(K)::Vector_3 &start, 
			 const typename CGAL_WRAP(K)::Vector_3 &end,
			 const K&)
{
    return  CGAL_NTS abs(wmult((K*)0, startwdist, end.hw())) -
            CGAL_NTS abs(wmult((K*)0, endwdist, start.hw()));
}


template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Segment_3 &seg1,
    const typename CGAL_WRAP(K)::Segment_3 &seg2,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 &start1 = seg1.start();
    const Point_3 &start2 = seg2.start();
    const Point_3 &end1 = seg1.end();
    const Point_3 &end2 = seg2.end();

    //    std::cout << "A" << std::endl;
    if (start1 == end1)
        return squared_distance(start1, seg2, k);
    if (start2 == end2)
        return squared_distance(start2, seg1, k);
    
    //    std::cout << "B" << std::endl;
    Vector_3 dir1, dir2, normal;
    dir1 = seg1.direction().vector();
    dir2 = seg2.direction().vector();
    normal = wcross(dir1, dir2, k);
    if (is_null(normal, k))
        return squared_distance_parallel(seg1, seg2, k);
    
    //    std::cout << "C" << std::endl;
    bool crossing1, crossing2;
    RT sdm_s1to2, sdm_e1to2, sdm_s2to1, sdm_e2to1;
    Vector_3 perpend1, perpend2, s2mins1, e2mins1, e1mins2;
    perpend1 = wcross(dir1, normal, k);
    perpend2 = wcross(dir2, normal, k);
    s2mins1 = start2-start1;
    e2mins1 = end2-start1;
    e1mins2 = end1-start2;
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
      //    std::cout << "crossing 1" << std::endl;
        if (crossing2) {
	  //	      std::cout << "crossing 2" << std::endl;
            return squared_distance_to_plane(normal, s2mins1, k);
        }
    
        RT dm;
        dm = _distance_measure_sub(
                  sdm_s2to1, sdm_e2to1, s2mins1, e2mins1, k);
        if (dm < RT(0)) {
	  //    std::cout << "X" << std::endl;
            return squared_distance(start2, seg1, k);
        } else {
            if (dm > RT(0)) {
	      //    std::cout << "Y" << std::endl;
                return squared_distance(end2, seg1, k);
            } else {
	      //    std::cout << "Z" << std::endl;
                // should not happen with exact arithmetic.
                return squared_distance_parallel(seg1, seg2, k);
            }
        }
    } else {
        if (crossing2) {
	  //    std::cout << "crossing 2" << std::endl;
            RT dm;
            dm =_distance_measure_sub(
                 sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
            if (dm < RT(0)) {
	      //	      std::cout << "X1" << std::endl;
                return squared_distance(start1, seg2, k);
            } else {
                if (dm > RT(0)) {
		  //    std::cout << "X2" << std::endl;
                    return squared_distance(end1, seg2, k);
                } else {
		  //    std::cout << "X3" << std::endl;
                    // should not happen with exact arithmetic.
                    return squared_distance_parallel(seg1, seg2, k);
                }
            }
        } else {
	  //    std::cout << "D" << std::endl;    
            FT min1, min2;
            RT dm;
            dm = _distance_measure_sub(
                     sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
            if (dm == RT(0)) // should not happen with exact arithmetic.
               return squared_distance_parallel(seg1, seg2, k);
            min1 = (dm < RT(0)) ?
                squared_distance(seg1.start(), seg2, k):
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
    const typename CGAL_WRAP(K)::Segment_3 &seg,
    const typename CGAL_WRAP(K)::Ray_3 &ray,
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
        if (!is_acute_angle(seg.start(), seg.end(), ray.start(), k))
            return squared_distance(seg.end(), ray.start(), k);
    } else {
        if (!is_acute_angle(seg.end(), seg.start(), ray.start(), k))
            return squared_distance(seg.start(), ray.start(), k);
    }
    return squared_distance(ray.start(), seg.supporting_line(), k);
}


template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Segment_3 &seg,
    const typename CGAL_WRAP(K)::Ray_3 &ray,
    const K& k)
{

  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 & ss = seg.start();
    const Point_3 & se = seg.end();
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
    ss_min_rs = ss-ray.start();
    se_min_rs = se-ray.start();
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
        return squared_distance(ray.start(), seg, k);
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
            min2 = squared_distance(ray.start(), seg, k);
            return (min1 < min2) ? min1 : min2;
        }
    }
}



template <class K>
inline
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Ray_3 & ray,
    const typename CGAL_WRAP(K)::Segment_3 & seg,
    const K& k)
{
    return squared_distance(seg, ray, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Segment_3 &seg,
    const typename CGAL_WRAP(K)::Line_3 &line,
    const K& k)
{

  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    const Point_3 &linepoint = line.point();
    const Point_3 &start = seg.start();
    const Point_3 &end = seg.end();

    if (start == end)
        return squared_distance(start, line, k);
    Vector_3 linedir = line.direction().vector();
    Vector_3 segdir = seg.direction().vector();
    Vector_3 normal = wcross(segdir, linedir, k);
    if (is_null(normal, k))
        return squared_distance_to_line(linedir, start-linepoint, k);

    bool crossing;
    RT sdm_ss2l, sdm_se2l;
    Vector_3 perpend2line, start_min_lp, end_min_lp;
    perpend2line = wcross(linedir, normal, k);
    start_min_lp = start-linepoint;
    end_min_lp = end-linepoint;
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
    const typename CGAL_WRAP(K)::Line_3 & line,
    const typename CGAL_WRAP(K)::Segment_3 & seg,
    const K& k)
{
    return squared_distance(seg, line, k);
}




template <class K>
typename K::FT
ray_ray_squared_distance_parallel(
    const typename CGAL_WRAP(K)::Vector_3 &ray1dir,
    const typename CGAL_WRAP(K)::Vector_3 &ray2dir,
    const typename CGAL_WRAP(K)::Vector_3 &s1_min_s2,
    const K& k)
{
    if (!is_acute_angle(ray2dir, s1_min_s2, k)) {
        bool same_direction;
        if (CGAL_NTS abs(ray1dir.hx()) > CGAL_NTS abs(ray1dir.hy())) {
            if (CGAL_NTS abs(ray1dir.hx()) > CGAL_NTS abs(ray1dir.hz()))
                same_direction =
                   (CGAL_NTS sign(ray1dir.hx()) == CGAL_NTS sign(ray2dir.hx()));
            else
                same_direction =
                   (CGAL_NTS sign(ray1dir.hz()) == CGAL_NTS sign(ray2dir.hz()));
        } else {
            if (CGAL_NTS abs(ray1dir.hy()) > CGAL_NTS abs(ray1dir.hz()))
                same_direction =
                   (CGAL_NTS sign(ray1dir.hy()) == CGAL_NTS sign(ray2dir.hy()));
            else
                same_direction =
                   (CGAL_NTS sign(ray1dir.hz()) == CGAL_NTS sign(ray2dir.hz()));
        }
        if (!same_direction)
            return (typename K::FT)(s1_min_s2*s1_min_s2);
    }
    return squared_distance_to_line(ray1dir, s1_min_s2, k);
}


template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Ray_3 &ray1,
    const typename CGAL_WRAP(K)::Ray_3 &ray2,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    const Point_3 & s1 = ray1.start();
    const Point_3 & s2 = ray2.start();
    Vector_3 dir1, dir2, normal;
    dir1 = ray1.direction().vector();
    dir2 = ray2.direction().vector();
    normal = wcross(dir1, dir2, k);
    Vector_3 s1_min_s2 = s1-s2;
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
        return squared_distance(ray2.start(), ray1, k);
    } else {
        if (crossing2) {
            return squared_distance(ray1.start(), ray2, k);
        } else {
          FT min1, min2;
            min1 = squared_distance(ray1.start(), ray2, k);
            min2 = squared_distance(ray2.start(), ray1, k);
            return (min1 < min2) ? min1 : min2;
        }
    }
}





template <class K>
extern typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Line_3 &line,
    const typename CGAL_WRAP(K)::Ray_3 &ray,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
    typedef typename K::RT RT;
    const Point_3 & rs =ray.start();
    Vector_3 raydir, linedir, normal;
    linedir = line.direction().vector();
    raydir = ray.direction().vector();
    normal = wcross(raydir, linedir, k);
    Vector_3 rs_min_lp = rs-line.point();
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
    const typename CGAL_WRAP(K)::Ray_3 & ray,
    const typename CGAL_WRAP(K)::Line_3 & line,
    const K& k)
{
    return squared_distance(line, ray, k);
}




template <class K>
typename K::FT
squared_distance(
    const typename CGAL_WRAP(K)::Line_3 &line1,
    const typename CGAL_WRAP(K)::Line_3 &line2,
    const K& k)
{
  typedef typename K::Vector_3 Vector_3;

    Vector_3 dir1, dir2, normal, diff;
    dir1 = line1.direction().vector();
    dir2 = line2.direction().vector();
    normal = wcross(dir1, dir2, k);
    diff = line2.point() - line1.point();
    if (is_null(normal, k))
        return squared_distance_to_line(dir2, diff, k);
    return squared_distance_to_plane(normal, diff, k);
}



} // namespace CGALi



template <class K>
inline
typename K::FT
squared_distance(const Point_3<K> &pt,
		 const Line_3<K> &line)
{
  return CGALi::squared_distance(pt, line, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> & line,
    const Point_3<K> & pt)
{
  return CGALi::squared_distance(pt, line, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Point_3<K> &pt,
    const Ray_3<K> &ray)
{
    return CGALi::squared_distance(pt, ray, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Point_3<K> & pt)
{
    return CGALi::squared_distance(pt, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Point_3<K> &pt,
    const Segment_3<K> &seg)
{
  return CGALi::squared_distance(pt, seg, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Segment_3<K> & seg,
    const Point_3<K> & pt)
{
    return CGALi::squared_distance(pt, seg, K());
}




template <class K>
inline
typename K::FT
squared_distance_parallel(
    const Segment_3<K> &seg1,
    const Segment_3<K> &seg2)
{
  return CGALi::squared_distance_parallel(seg1, seg2, K());
}




template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K> &seg1,
		 const Segment_3<K> &seg2)
{
  return CGALi::squared_distance(seg1, seg2, K());
}






template <class K>
inline
typename K::FT
squared_distance_parallel(
    const Segment_3<K> &seg,
    const Ray_3<K> &ray)
{
  return CGALi::squared_distance_parallel(ray,seg, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Segment_3<K> &seg,
    const Ray_3<K> &ray)
{
  return CGALi::squared_distance(ray, seg, K());
}



template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Segment_3<K> & seg)
{
    return squared_distance(seg, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Segment_3<K> &seg,
    const Line_3<K> &line)
{
  return CGALi::squared_distance(seg, line, K());
}


template <class K>
inline 
typename K::FT
squared_distance(
    const Line_3<K> & line,
    const Segment_3<K> & seg)
{
    return CGALi::squared_distance(seg, line, K());
}




template <class K>
inline
typename K::FT
ray_ray_squared_distance_parallel(
    const Vector_3<K> &ray1dir,
    const Vector_3<K> &ray2dir,
    const Vector_3<K> &s1_min_s2)
{
 
  return CGALi::ray_ray_squared_distance_parallel(ray1dir, ray2dir, s1_min_s2, K());
}

template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> &ray1,
    const Ray_3<K> &ray2)
{
  return CGALi::squared_distance(ray1, ray2, K());
}





template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> &line,
    const Ray_3<K> &ray)
{
  return CGALi::squared_distance(line, ray, K());
}


template <class K>
inline
typename K::FT
squared_distance(
    const Ray_3<K> & ray,
    const Line_3<K> & line)
{
    return CGALi::squared_distance(line, ray, K());
}




template <class K>
inline
typename K::FT
squared_distance(
    const Line_3<K> &line1,
    const Line_3<K> &line2)
{  
    return CGALi::squared_distance(line1, line2, K());
}



CGAL_END_NAMESPACE


#endif
