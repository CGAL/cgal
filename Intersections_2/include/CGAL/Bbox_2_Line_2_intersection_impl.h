// Copyright (c) 2000
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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> Lcart;


namespace CGAL {

class Bbox_2_Line_2_pair_impl
{
public:
    Bbox_2_Line_2_pair_impl() {}
    Bbox_2_Line_2_pair_impl(Bbox_2 const &bb, Lcart::Line_2 const &line)
        : _bbox(bb), _line(line), _known(false) {}
    Bbox_2 _bbox;
    Lcart::Line_2 _line;
    mutable bool                     _known;
    mutable Bbox_2_Line_2_pair::Intersection_results     _result;
    mutable double                   _min, _max;
};

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair::~Bbox_2_Line_2_pair()
{
    delete pimpl;
}

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair::Bbox_2_Line_2_pair()
{
    pimpl = new Bbox_2_Line_2_pair_impl;
    pimpl->_known = false;
}

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair::Bbox_2_Line_2_pair(Bbox_2_Line_2_pair const &o)
{
    pimpl = new Bbox_2_Line_2_pair_impl(*o.pimpl);
}

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair::Bbox_2_Line_2_pair(
    Bbox_2 const &bbox, double line_a, double line_b, double line_c)
{
    pimpl = new Bbox_2_Line_2_pair_impl(bbox,
              Lcart::Line_2(line_a, line_b, line_c));
}

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair &
Bbox_2_Line_2_pair::operator=(Bbox_2_Line_2_pair const &o)
{
    *pimpl = *o.pimpl;
    return *this;
}

CGAL_INLINE_FUNCTION
Bbox_2_Line_2_pair::Intersection_results
Bbox_2_Line_2_pair::intersection_type() const
{
    if (pimpl->_known)
        return pimpl->_result;
    // The non const this pointer is used to cast away const.
    pimpl->_known = true;
    const Lcart::Point_2 &ref_point = pimpl->_line.point();
    const Lcart::Vector_2 &dir =
                               pimpl->_line.direction().to_vector();
    bool to_infinity = true;
// first on x value
    if (dir.x() == 0.0) {
        if (ref_point.x() < pimpl->_bbox.xmin()) {
            pimpl->_result = NO_INTERSECTION;
            return pimpl->_result;
        }
        if (ref_point.x() > pimpl->_bbox.xmax()) {
            pimpl->_result = NO_INTERSECTION;
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
                pimpl->_result = NO_INTERSECTION;
                return pimpl->_result;
            }
        }
        to_infinity = false;
    }
// now on y value
    if (dir.y() == 0.0) {
        if (ref_point.y() < pimpl->_bbox.ymin()) {
            pimpl->_result = NO_INTERSECTION;
            return pimpl->_result;
        }
        if (ref_point.y() > pimpl->_bbox.ymax()) {
            pimpl->_result = NO_INTERSECTION;
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
                pimpl->_result = NO_INTERSECTION;
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

CGAL_INLINE_FUNCTION
bool
Bbox_2_Line_2_pair::intersection(
    double &x1, double &y1, double &x2, double &y2) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != SEGMENT)
        return false;
    Lcart::Point_2 p1(pimpl->_line.point()
                + pimpl->_min*pimpl->_line.direction().to_vector());
    Lcart::Point_2 p2(pimpl->_line.point()
                + pimpl->_max*pimpl->_line.direction().to_vector());
    x1 = p1.x();
    y1 = p1.y();
    x2 = p2.x();
    y2 = p2.y();
    return true;
}

CGAL_INLINE_FUNCTION
bool
Bbox_2_Line_2_pair::intersection(
    double &x, double &y) const
{
    if (!pimpl->_known)
        intersection_type();
    if (pimpl->_result != POINT)
        return false;
    Lcart::Point_2 pt(pimpl->_line.point()
        + pimpl->_min*pimpl->_line.direction().to_vector());
    x = pt.x();
    y = pt.y();
    return true;
}

} //namespace CGAL
