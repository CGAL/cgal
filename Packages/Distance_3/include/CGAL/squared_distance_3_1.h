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



template <class R>
typename R::FT
squared_distance(
    const Point_3<R> &pt,
    const Line_3<R> &line)
{
    Vector_3<R> dir(line.direction().vector());
    Vector_3<R> diff = pt - line.point();
    return squared_distance_to_line(dir, diff);
}


template <class R>
inline typename R::FT
squared_distance(
    const Line_3<R> & line,
    const Point_3<R> & pt)
{
    return squared_distance(pt, line);
}


template <class R>
extern typename R::FT
squared_distance(
    const Point_3<R> &pt,
    const Ray_3<R> &ray)
{
    Vector_3<R> diff = pt-ray.start();
    const Vector_3<R> &dir = ray.direction().vector();
    if (!is_acute_angle(dir,diff) )
        return (typename R::FT)(diff*diff);
    return squared_distance_to_line(dir, diff);
}


template <class R>
inline typename R::FT
squared_distance(
    const Ray_3<R> & ray,
    const Point_3<R> & pt)
{
    return squared_distance(pt, ray);
}


template <class R>
typename R::FT
squared_distance(
    const Point_3<R> &pt,
    const Segment_3<R> &seg)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    // assert that the segment is valid (non zero length).
    Vector_3<R> diff = pt-seg.start();
    Vector_3<R> segvec = seg.end()-seg.start();
    RT d = wdot(diff,segvec);
    if (d <= (RT)0)
        return (FT(diff*diff));
    RT e = wdot(segvec,segvec);
    if (wmult((R*)0 ,d, segvec.hw()) > wmult((R*)0, e, diff.hw()))
        return squared_distance(pt, seg.end());
    return squared_distance_to_line(segvec, diff);
}


template <class R>
inline typename R::FT
squared_distance(
    const Segment_3<R> & seg,
    const Point_3<R> & pt)
{
    return squared_distance(pt, seg);
}




template <class R>
typename R::FT
squared_distance_parallel(
    const Segment_3<R> &seg1,
    const Segment_3<R> &seg2)
{
    bool same_direction;
    const Vector_3<R> &dir1 = seg1.direction().vector();
    const Vector_3<R> &dir2 = seg2.direction().vector();
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
        if (!is_acute_angle(seg1.start(), seg1.end(), seg2.start()))
            return squared_distance(seg1.end(), seg2.start());
        if (!is_acute_angle(seg1.end(), seg1.start(), seg2.end()))
            return squared_distance(seg1.start(), seg2.end());
    } else {
        if (!is_acute_angle(seg1.start(), seg1.end(), seg2.end()))
            return squared_distance(seg1.end(), seg2.end());
        if (!is_acute_angle(seg1.end(), seg1.start(), seg2.start()))
            return squared_distance(seg1.start(), seg2.start());
    }
    return squared_distance(seg2.start(), seg1.supporting_line());
}



template <class RT, class R>
RT _distance_measure_sub(RT startwdist, RT endwdist,
const Vector_3<R> &start, const Vector_3<R> &end
)
{
    return  CGAL_NTS abs(wmult((R*)0, startwdist, end.hw())) -
            CGAL_NTS abs(wmult((R*)0, endwdist, start.hw()));
}


template <class R>
typename R::FT
squared_distance(
    const Segment_3<R> &seg1,
    const Segment_3<R> &seg2)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    const Point_3<R> &start1 = seg1.start();
    const Point_3<R> &start2 = seg2.start();
    const Point_3<R> &end1 = seg1.end();
    const Point_3<R> &end2 = seg2.end();

    if (start1 == end1)
        return squared_distance(start1, seg2);
    if (start2 == end2)
        return squared_distance(start2, seg1);
    
    Vector_3<R> dir1, dir2, normal;
    dir1 = seg1.direction().vector();
    dir2 = seg2.direction().vector();
    normal = wcross(dir1, dir2);
    if (is_null(normal))
        return squared_distance_parallel(seg1, seg2);
    
    bool crossing1, crossing2;
    RT sdm_s1to2, sdm_e1to2, sdm_s2to1, sdm_e2to1;
    Vector_3<R> perpend1, perpend2, s2mins1, e2mins1, e1mins2;
    perpend1 = wcross(dir1, normal);
    perpend2 = wcross(dir2, normal);
    s2mins1 = start2-start1;
    e2mins1 = end2-start1;
    e1mins2 = end1-start2;
    sdm_s1to2 = -RT(wdot(perpend2, s2mins1));
    sdm_e1to2 = wdot(perpend2, e1mins2);
    sdm_s2to1 = wdot(perpend1, s2mins1);
    sdm_e2to1 = wdot(perpend1, e2mins1);
    
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
            return squared_distance_to_plane(normal, s2mins1);
        }
    
        RT dm;
        dm = _distance_measure_sub(
                  sdm_s2to1, sdm_e2to1, s2mins1, e2mins1);
        if (dm < RT(0)) {
            return squared_distance(start2, seg1);
        } else {
            if (dm > RT(0)) {
                return squared_distance(end2, seg1);
            } else {
                // should not happen with exact arithmetic.
                return squared_distance_parallel(seg1, seg2);
            }
        }
    } else {
        if (crossing2) {
            RT dm;
            dm =_distance_measure_sub(
                 sdm_s1to2, sdm_e1to2, s2mins1, e1mins2);
            if (dm < RT(0)) {
                return squared_distance(start1, seg2);
            } else {
                if (dm > RT(0)) {
                    return squared_distance(end1, seg2);
                } else {
                    // should not happen with exact arithmetic.
                    return squared_distance_parallel(seg1, seg2);
                }
            }
        } else {
    
            FT min1, min2;
            RT dm;
            dm = _distance_measure_sub(
                     sdm_s1to2, sdm_e1to2, s2mins1, e1mins2);
            if (dm == RT(0)) // should not happen with exact arithmetic.
               return squared_distance_parallel(seg1, seg2);
            min1 = (dm < RT(0)) ?
                squared_distance(seg1.start(), seg2):
                squared_distance(end1, seg2);
            dm = _distance_measure_sub(
                     sdm_s2to1, sdm_e2to1, s2mins1, e2mins1);
            if (dm == RT(0)) // should not happen with exact arithmetic.
                return squared_distance_parallel(seg1, seg2);
            min2 = (dm < RT(0)) ?
                squared_distance(start2, seg1):
                squared_distance(end2, seg1);
            return (min1 < min2) ? min1 : min2;
        }
    }
    
}






template <class R>
typename R::FT
squared_distance_parallel(
    const Segment_3<R> &seg,
    const Ray_3<R> &ray)
{
    bool same_direction;
    const Vector_3<R> &dir1 = seg.direction().vector();
    const Vector_3<R> &dir2 = ray.direction().vector();
    if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
        same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
    } else {
        same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
    }
    if (same_direction) {
        if (!is_acute_angle(seg.start(), seg.end(), ray.start()))
            return squared_distance(seg.end(), ray.start());
    } else {
        if (!is_acute_angle(seg.end(), seg.start(), ray.start()))
            return squared_distance(seg.start(), ray.start());
    }
    return squared_distance(ray.start(), seg.supporting_line());
}


template <class R>
typename R::FT
squared_distance(
    const Segment_3<R> &seg,
    const Ray_3<R> &ray)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    const Point_3<R> & ss = seg.start();
    const Point_3<R> & se = seg.end();
    if (ss == se)
        return squared_distance(ss, ray);
    Vector_3<R> raydir, segdir, normal;
    raydir = ray.direction().vector();
    segdir = seg.direction().vector();
    normal = wcross(segdir, raydir);
    if (is_null(normal))
        return squared_distance_parallel(seg, ray);

    bool crossing1, crossing2;
    RT sdm_ss2r, sdm_se2r, sdm_rs2s, sdm_re2s;
    Vector_3<R> perpend2seg, perpend2ray, ss_min_rs, se_min_rs;
    perpend2seg = wcross(segdir, normal);
    perpend2ray = wcross(raydir, normal);
    ss_min_rs = ss-ray.start();
    se_min_rs = se-ray.start();
    sdm_ss2r = wdot(perpend2ray, ss_min_rs);
    sdm_se2r = wdot(perpend2ray, se_min_rs);
    if (sdm_ss2r < RT(0)) {
        crossing1 = (sdm_se2r >= RT(0));
    } else {
        if (sdm_se2r <= RT(0)) {
            crossing1 = true;
        } else {
            crossing1 = (sdm_ss2r == RT(0));
        }
    }

    sdm_rs2s = -RT(wdot(perpend2seg, ss_min_rs));
    sdm_re2s = wdot(perpend2seg, raydir);
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
            return squared_distance_to_plane(normal, ss_min_rs);
        }
        return squared_distance(ray.start(), seg);
    } else {
        if (crossing2) {
            RT dm;
            dm = _distance_measure_sub(
                    sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs);
            if (dm < RT(0)) {
                return squared_distance(ss, ray);
            } else {
                if (dm > RT(0)) {
                    return squared_distance(se, ray);
                } else {
                    // parallel, should not happen (no crossing)
                    return squared_distance_parallel(seg, ray);
                }
            }
        } else {
            FT min1, min2;
            RT dm;
            dm = _distance_measure_sub(
                    sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs);
            if (dm == RT(0))
                return squared_distance_parallel(seg, ray);
            min1 = (dm < RT(0))
                 ? squared_distance(ss, ray)
                 : squared_distance(se, ray);
            min2 = squared_distance(ray.start(), seg);
            return (min1 < min2) ? min1 : min2;
        }
    }
}



template <class R>
inline typename R::FT
squared_distance(
    const Ray_3<R> & ray,
    const Segment_3<R> & seg)
{
    return squared_distance(seg, ray);
}


template <class R>
typename R::FT
squared_distance(
    const Segment_3<R> &seg,
    const Line_3<R> &line)
{
    typedef typename R::RT RT;
    const Point_3<R> &linepoint = line.point();
    const Point_3<R> &start = seg.start();
    const Point_3<R> &end = seg.end();

    if (start == end)
        return squared_distance(start, line);
    Vector_3<R> linedir = line.direction().vector();
    Vector_3<R> segdir = seg.direction().vector();
    Vector_3<R> normal = wcross(segdir, linedir);
    if (is_null(normal))
        return squared_distance_to_line(linedir, start-linepoint);

    bool crossing;
    RT sdm_ss2l, sdm_se2l;
    Vector_3<R> perpend2line, start_min_lp, end_min_lp;
    perpend2line = wcross(linedir, normal);
    start_min_lp = start-linepoint;
    end_min_lp = end-linepoint;
    sdm_ss2l = wdot(perpend2line, start_min_lp);
    sdm_se2l = wdot(perpend2line, end_min_lp);
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
        return squared_distance_to_plane(normal, start_min_lp);
    } else {
        RT dm;
        dm = _distance_measure_sub(
                sdm_ss2l, sdm_se2l, start_min_lp, end_min_lp);
        if (dm <= RT(0)) {
            return squared_distance_to_line(linedir, start_min_lp);
        } else {
            return squared_distance_to_line(linedir, end_min_lp);
        }
    }
}


template <class R>
inline typename R::FT
squared_distance(
    const Line_3<R> & line,
    const Segment_3<R> & seg)
{
    return squared_distance(seg, line);
}




template <class R>
typename R::FT
ray_ray_squared_distance_parallel(
    const Vector_3<R> &ray1dir,
    const Vector_3<R> &ray2dir,
    const Vector_3<R> &s1_min_s2)
{
    if (!is_acute_angle(ray2dir, s1_min_s2)) {
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
            return (typename R::FT)(s1_min_s2*s1_min_s2);
    }
    return squared_distance_to_line(ray1dir, s1_min_s2);
}


template <class R>
typename R::FT
squared_distance(
    const Ray_3<R> &ray1,
    const Ray_3<R> &ray2)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    const Point_3<R> & s1 = ray1.start();
    const Point_3<R> & s2 = ray2.start();
    Vector_3<R> dir1, dir2, normal;
    dir1 = ray1.direction().vector();
    dir2 = ray2.direction().vector();
    normal = wcross(dir1, dir2);
    Vector_3<R> s1_min_s2 = s1-s2;
    if (is_null(normal))
        return ray_ray_squared_distance_parallel(dir1, dir2, s1_min_s2);

    bool crossing1, crossing2;
    RT sdm_s1_2, sdm_s2_1;
    Vector_3<R> perpend1, perpend2;
    perpend1 = wcross(dir1, normal);
    perpend2 = wcross(dir2, normal);

    sdm_s1_2 = wdot(perpend2, s1_min_s2);
    if (sdm_s1_2 < RT(0)) {
        crossing1 = (RT(wdot(perpend2, dir1)) >= RT(0));
    } else {
        if (RT(wdot(perpend2, dir1)) <= RT(0)) {
            crossing1 = true;
        } else {
            crossing1 = (sdm_s1_2 == RT(0));
        }
    }
    sdm_s2_1 = -RT(wdot(perpend1, s1_min_s2));
    if (sdm_s2_1 < RT(0)) {
        crossing2 = (RT(wdot(perpend1, dir2)) >= RT(0));
    } else {
        if (RT(wdot(perpend1, dir2)) <= RT(0)) {
            crossing2 = true;
        } else {
            crossing2 = (sdm_s2_1 == RT(0));
        }
    }
    if (crossing1) {
        if (crossing2)
            return squared_distance_to_plane(normal, s1_min_s2);
        return squared_distance(ray2.start(), ray1);
    } else {
        if (crossing2) {
            return squared_distance(ray1.start(), ray2);
        } else {
          FT min1, min2;
            min1 = squared_distance(ray1.start(), ray2);
            min2 = squared_distance(ray2.start(), ray1);
            return (min1 < min2) ? min1 : min2;
        }
    }
}





template <class R>
extern typename R::FT
squared_distance(
    const Line_3<R> &line,
    const Ray_3<R> &ray)
{
    typedef typename R::RT RT;
    const Point_3<R> & rs =ray.start();
    Vector_3<R> raydir, linedir, normal;
    linedir = line.direction().vector();
    raydir = ray.direction().vector();
    normal = wcross(raydir, linedir);
    Vector_3<R> rs_min_lp = rs-line.point();
    if (is_null(normal))
        return squared_distance_to_line(linedir, rs_min_lp);

    bool crossing;
    RT sdm_sr_l;
    Vector_3<R> perpend2l;
    perpend2l = wcross(linedir, normal);

    sdm_sr_l = wdot(perpend2l, rs_min_lp);
    if (sdm_sr_l < RT(0)) {
        crossing = (RT(wdot(perpend2l, raydir)) >= RT(0));
    } else {
        if (RT(wdot(perpend2l, raydir)) <= RT(0)) {
            crossing = true;
        } else {
            crossing = (sdm_sr_l == RT(0));
        }
    }

    if (crossing)
        return squared_distance_to_plane(normal, rs_min_lp);
    else
        return squared_distance_to_line(linedir, rs_min_lp);
}


template <class R>
inline typename R::FT
squared_distance(
    const Ray_3<R> & ray,
    const Line_3<R> & line)
{
    return squared_distance(line, ray);
}




template <class R>
typename R::FT
squared_distance(
    const Line_3<R> &line1,
    const Line_3<R> &line2)
{
    Vector_3<R> dir1, dir2, normal, diff;
    dir1 = line1.direction().vector();
    dir2 = line2.direction().vector();
    normal = wcross(dir1, dir2);
    diff = line2.point() - line1.point();
    if (is_null(normal))
        return squared_distance_to_line(dir2, diff);
    return squared_distance_to_plane(normal, diff);
}



CGAL_END_NAMESPACE


#endif
