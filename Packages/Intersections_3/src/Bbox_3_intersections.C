// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : src/Bbox_3_intersections.C
// source        : web/intersection_3.fw
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Saarbruecken
//
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/bbox_intersection_3.h>


CGAL_BEGIN_NAMESPACE

Object
intersection_bl(const Bbox_3 &box,
        double lpx, double lpy, double lpz,
        double ldx, double ldy, double ldz,
        bool min_infinite, bool max_infinite)
{
    typedef Cartesian<double> R_cd;
    Point_3<R_cd> ref_point(lpx,lpy, lpz);
    Vector_3<R_cd> dir(ldx, ldy, ldz);
    double seg_min = 0.0, seg_max = 1.0;
// first on x value
    if (dir.x() == 0.0) {
        if (ref_point.x() < box.xmin())
            return Object();
        if (ref_point.x() > box.xmax())
            return Object();
    } else {
        double newmin, newmax;
        if (dir.x() > 0.0) {
            newmin = (box.xmin()-ref_point.x())/dir.x();
            newmax = (box.xmax()-ref_point.x())/dir.x();
        } else {
            newmin = (box.xmax()-ref_point.x())/dir.x();
            newmax = (box.xmin()-ref_point.x())/dir.x();
        }
        if (min_infinite) {
            min_infinite = false;
            seg_min = newmin;
        } else {
            if (newmin > seg_min)
                seg_min = newmin;
        }
        if (max_infinite) {
            max_infinite = false;
            seg_max = newmax;
        } else {
            if (newmax < seg_max)
                seg_max = newmax;
        }
        if (seg_max < seg_min)
            return Object();
    }
// now on y value
    if (dir.y() == 0.0) {
        if (ref_point.y() < box.ymin())
            return Object();
        if (ref_point.y() > box.ymax())
            return Object();
    } else {
        double newmin, newmax;
        if (dir.y() > 0.0) {
            newmin = (box.ymin()-ref_point.y())/dir.y();
            newmax = (box.ymax()-ref_point.y())/dir.y();
        } else {
            newmin = (box.ymax()-ref_point.y())/dir.y();
            newmax = (box.ymin()-ref_point.y())/dir.y();
        }
        if (min_infinite) {
            min_infinite = false;
            seg_min = newmin;
        } else {
            if (newmin > seg_min)
                seg_min = newmin;
        }
        if (max_infinite) {
            max_infinite = false;
            seg_max = newmax;
        } else {
            if (newmax < seg_max)
                seg_max = newmax;
        }
        if (seg_max < seg_min)
            return Object();
    }
// now on z value
    if (dir.z() == 0.0) {
        if (ref_point.z() < box.zmin())
            return Object();
        if (ref_point.z() > box.zmax())
            return Object();
    } else {
        double newmin, newmax;
        if (dir.z() > 0.0) {
            newmin = (box.zmin()-ref_point.z())/dir.z();
            newmax = (box.zmax()-ref_point.z())/dir.z();
        } else {
            newmin = (box.zmax()-ref_point.z())/dir.z();
            newmax = (box.zmin()-ref_point.z())/dir.z();
        }
        if (min_infinite) {
            min_infinite = false;
            seg_min = newmin;
        } else {
            if (newmin > seg_min)
                seg_min = newmin;
        }
        if (max_infinite) {
            max_infinite = false;
            seg_max = newmax;
        } else {
            if (newmax < seg_max)
                seg_max = newmax;
        }
        if (seg_max < seg_min)
            return Object();
    }
    if (min_infinite || max_infinite) {
        seg_max = 0.0;
        CGAL_kernel_assertion_msg(true,
            "Zero direction vector of line detected.");
    }
    if (seg_max == seg_min)
        return make_object(ref_point + dir*seg_max);
    return make_object(Segment_3<R_cd>(
            ref_point + dir*seg_min, ref_point + dir*seg_max));
}

CGAL_END_NAMESPACE

