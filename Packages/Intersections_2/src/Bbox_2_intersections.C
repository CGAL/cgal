// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : src/Bbox_2_intersections.C
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bbox_2_Line_2_intersection.h>
#include <CGAL/Ray_2_Bbox_2_intersection.h>
typedef CGAL::Simple_cartesian<double> Rcart;


CGAL_BEGIN_NAMESPACE

class Bbox_2_Line_2_pair_impl
{
public:
    Bbox_2_Line_2_pair_impl() {}
    Bbox_2_Line_2_pair_impl(Bbox_2 const &bb, Rcart::Line_2 const &line)
        : _bbox(bb), _line(line), _known(false) {}
    Bbox_2 _bbox;
    Rcart::Line_2 _line;
    mutable bool                     _known;
    mutable Bbox_2_Line_2_pair::Intersection_results     _result;
    mutable double                   _min, _max;
};

Bbox_2_Line_2_pair::~Bbox_2_Line_2_pair()
{
    delete pimpl;
}

Bbox_2_Line_2_pair::Bbox_2_Line_2_pair()
{
    pimpl = new Bbox_2_Line_2_pair_impl;
    pimpl->_known = false;
}

Bbox_2_Line_2_pair::Bbox_2_Line_2_pair(Bbox_2_Line_2_pair const &o)
{
    pimpl = new Bbox_2_Line_2_pair_impl(*o.pimpl);
}

Bbox_2_Line_2_pair::Bbox_2_Line_2_pair(
    Bbox_2 const &bbox, double line_a, double line_b, double line_c)
{
    pimpl = new Bbox_2_Line_2_pair_impl(bbox,
              Rcart::Line_2(line_a, line_b, line_c));
}

Bbox_2_Line_2_pair &
Bbox_2_Line_2_pair::operator=(Bbox_2_Line_2_pair const &o)
{
    *pimpl = *o.pimpl;
    return *this;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
Bbox_2_Line_2_pair::Intersection_results
Bbox_2_Line_2_pair::intersection_type() const
{
    if (pimpl->_known)
        return pimpl->_result;
    // The non const this pointer is used to cast away const.
    pimpl->_known = true;
    const Rcart::Point_2 &ref_point = pimpl->_line.point();
    const Rcart::Vector_2 &dir =
                               pimpl->_line.direction().to_vector();
    bool to_infinity = true;
// first on x value
    if (dir.x() == 0.0) {
        if (ref_point.x() < pimpl->_bbox.xmin()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
        if (ref_point.x() > pimpl->_bbox.xmax()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
    } else {
        double newmin, newmax;
        if (dir.x() > 0.0) {
            newmin = (pimpl->_bbox.xmin()-ref_point.x())/dir.x();
            newmax = (pimpl->_bbox.xmax()-ref_point.x())/dir.x();
        } else {
            newmin = (pimpl->_bbox.xmax()-ref_point.x())/dir.x();
            newmax = (pimpl->_bbox.xmin()-ref_point.x())/dir.x();
        }
        if (to_infinity) {
            pimpl->_min = newmin;
            pimpl->_max = newmax;
        } else {
            if (newmin > pimpl->_min)
                pimpl->_min = newmin;
            if (newmax < pimpl->_max)
                pimpl->_max = newmax;
            if (pimpl->_max < pimpl->_min) {
                pimpl->_result = NO;
                return pimpl->_result;
            }
        }
        to_infinity = false;
    }
// now on y value
    if (dir.y() == 0.0) {
        if (ref_point.y() < pimpl->_bbox.ymin()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
        if (ref_point.y() > pimpl->_bbox.ymax()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
    } else {
        double newmin, newmax;
        if (dir.y() > 0.0) {
            newmin = (pimpl->_bbox.ymin()-ref_point.y())/dir.y();
            newmax = (pimpl->_bbox.ymax()-ref_point.y())/dir.y();
        } else {
            newmin = (pimpl->_bbox.ymax()-ref_point.y())/dir.y();
            newmax = (pimpl->_bbox.ymin()-ref_point.y())/dir.y();
        }
        if (to_infinity) {
            pimpl->_min = newmin;
            pimpl->_max = newmax;
        } else {
            if (newmin > pimpl->_min)
                pimpl->_min = newmin;
            if (newmax < pimpl->_max)
                pimpl->_max = newmax;
            if (pimpl->_max < pimpl->_min) {
                pimpl->_result = NO;
                return pimpl->_result;
            }
        }
        to_infinity = false;
    }
    CGAL_kernel_assertion(!to_infinity);
    if (pimpl->_max == pimpl->_min) {
        pimpl->_result = POINT;
        return pimpl->_result;
    }
    pimpl->_result = SEGMENT;
    return pimpl->_result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

bool
Bbox_2_Line_2_pair::intersection(
    double &x1, double &y1, double &x2, double &y2) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != SEGMENT)
        return false;
    Rcart::Point_2 p1(pimpl->_line.point()
                + pimpl->_min*pimpl->_line.direction().to_vector());
    Rcart::Point_2 p2(pimpl->_line.point()
                + pimpl->_max*pimpl->_line.direction().to_vector());
    x1 = p1.x();
    y1 = p1.y();
    x2 = p2.x();
    y2 = p2.y();
    return true;
}

bool
Bbox_2_Line_2_pair::intersection(
    double &x, double &y) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != POINT)
        return false;
    Rcart::Point_2 pt(pimpl->_line.point()
        + pimpl->_min*pimpl->_line.direction().to_vector());
    x = pt.x();
    y = pt.y();
    return true;
}

CGAL_END_NAMESPACE






CGAL_BEGIN_NAMESPACE

class Bbox_2_Ray_2_pair_impl
{
public:
    Bbox_2_Ray_2_pair_impl():_known(false) {}
    Bbox_2_Ray_2_pair_impl(Bbox_2 const &bbox, Rcart::Point_2 const &pt,
                Rcart::Vector_2 const &dir)
        :_box(bbox), _known(false), _ref_point(pt), _dir(dir), _min(0.0) {}
    Ray_2< Rcart > _ray;
    Bbox_2 _box;
    bool _known;
    Bbox_2_Ray_2_pair::Intersection_results _result;
    Rcart::Point_2 _ref_point;
    Rcart::Vector_2 _dir;
    double _min, _max;
};

Bbox_2_Ray_2_pair::~Bbox_2_Ray_2_pair()
{
    delete pimpl;
}

Bbox_2_Ray_2_pair::Bbox_2_Ray_2_pair()
{
    pimpl = new Bbox_2_Ray_2_pair_impl;
}

Bbox_2_Ray_2_pair::Bbox_2_Ray_2_pair(Bbox_2_Ray_2_pair const &o)
{
    pimpl = new Bbox_2_Ray_2_pair_impl(*o.pimpl);
}

Bbox_2_Ray_2_pair::Bbox_2_Ray_2_pair(
    Bbox_2 const &bbox, double x, double y, double dx, double dy)
{
    pimpl = new Bbox_2_Ray_2_pair_impl(bbox,
                    Rcart::Point_2(x,y), Rcart::Vector_2(dx,dy));
}

Bbox_2_Ray_2_pair &
Bbox_2_Ray_2_pair::operator=(Bbox_2_Ray_2_pair const &o)
{
    *pimpl = *o.pimpl;
    return *this;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

Bbox_2_Ray_2_pair::Intersection_results
Bbox_2_Ray_2_pair::intersection_type() const
{
    if (pimpl->_known)
        return pimpl->_result;
    pimpl->_known = true;
    bool to_infinity = true;
// first on x value
    if (pimpl->_dir.x() == 0.0) {
        if (pimpl->_ref_point.x() < pimpl->_box.xmin()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
        if (pimpl->_ref_point.x() > pimpl->_box.xmax()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
    } else {
        double newmin, newmax;
        if (pimpl->_dir.x() > 0.0) {
            newmin = (pimpl->_box.xmin()-pimpl->_ref_point.x())/pimpl->_dir.x();
            newmax = (pimpl->_box.xmax()-pimpl->_ref_point.x())/pimpl->_dir.x();
        } else {
            newmin = (pimpl->_box.xmax()-pimpl->_ref_point.x())/pimpl->_dir.x();
            newmax = (pimpl->_box.xmin()-pimpl->_ref_point.x())/pimpl->_dir.x();
        }
        if (newmin > pimpl->_min)
            pimpl->_min = newmin;
        if (to_infinity) {
            pimpl->_max = newmax;
        } else {
            if (newmax < pimpl->_max)
                pimpl->_max = newmax;
        }
        if (pimpl->_max < pimpl->_min){
            pimpl->_result = NO;
            return pimpl->_result;
        }
        to_infinity = false;
    }
// now on y value
    if (pimpl->_dir.y() == 0.0) {
        if (pimpl->_ref_point.y() < pimpl->_box.ymin()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
        if (pimpl->_ref_point.y() > pimpl->_box.ymax()) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
    } else {
        double newmin, newmax;
        if (pimpl->_dir.y() > 0.0) {
            newmin = (pimpl->_box.ymin()-pimpl->_ref_point.y())/pimpl->_dir.y();
            newmax = (pimpl->_box.ymax()-pimpl->_ref_point.y())/pimpl->_dir.y();
        } else {
            newmin = (pimpl->_box.ymax()-pimpl->_ref_point.y())/pimpl->_dir.y();
            newmax = (pimpl->_box.ymin()-pimpl->_ref_point.y())/pimpl->_dir.y();
        }
        if (newmin > pimpl->_min)
            pimpl->_min = newmin;
        if (to_infinity) {
            pimpl->_max = newmax;
        } else {
            if (newmax < pimpl->_max)
                pimpl->_max = newmax;
        }
        if (pimpl->_max < pimpl->_min) {
            pimpl->_result = NO;
            return pimpl->_result;
        }
        to_infinity = false;
    }
    CGAL_kernel_assertion(!to_infinity);
    if (pimpl->_max == pimpl->_min) {
        pimpl->_result = POINT;
        return pimpl->_result;
    }
    pimpl->_result = SEGMENT;
    return pimpl->_result;
}


bool Bbox_2_Ray_2_pair::
intersection(double &x1, double &y1, double &x2, double &y2) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != SEGMENT)
        return false;
    Rcart::Point_2 p1(pimpl->_ref_point + pimpl->_min*pimpl->_dir);
    Rcart::Point_2 p2(pimpl->_ref_point + pimpl->_max*pimpl->_dir);
    x1 = p1.x();
    y1 = p1.y();
    x2 = p2.x();
    y2 = p2.y();
    return true;
}

bool Bbox_2_Ray_2_pair::intersection(double &x, double &y) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != POINT)
        return false;
    Rcart::Point_2 pt = pimpl->_ref_point + pimpl->_min*pimpl->_dir;
    x = pt.x();
    y = pt.y();
    return true;
}


bool do_intersect_ray_2(
    const Bbox_2 &box, double x, double y, double dx, double dy)
{
    Bbox_2_Ray_2_pair pair(box, x, y, dx, dy);
    return pair.intersection_type() != Bbox_2_Ray_2_pair::NO;
}

CGAL_END_NAMESPACE



