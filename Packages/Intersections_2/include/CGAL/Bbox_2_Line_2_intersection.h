
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
// file          : include/CGAL/Bbox_2_Line_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken
//
// ============================================================================


#ifndef CGAL_BBOX_2_LINE_2_INTERSECTION_H
#define CGAL_BBOX_2_LINE_2_INTERSECTION_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Line_2.h>
//#include <CGAL/Segment_2.h>
//#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

class Bbox_2_Line_2_pair_impl;

class Bbox_2_Line_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Bbox_2_Line_2_pair() ;
    Bbox_2_Line_2_pair(Bbox_2_Line_2_pair const &);
    Bbox_2_Line_2_pair(Bbox_2 const &bbox,
                            double line_a, double line_b, double line_c);
    ~Bbox_2_Line_2_pair() ;
    Bbox_2_Line_2_pair &operator=(Bbox_2_Line_2_pair const &o);
    // set_bbox(Bbox_2 const &bbox);
    // set_line(double line_a, double line_b, double line_c);
    Intersection_results intersection_type() const;
    bool intersection(double &x, double &y) const;
    bool intersection(double &x1, double &y1, double &x2, double &y2) const;
protected:
    Bbox_2_Line_2_pair_impl *pimpl;
};

template <class Line>
Bbox_2_Line_2_pair intersection_computer_line_2(
    Bbox_2 const &bbox, Line const &line)
{
    return Bbox_2_Line_2_pair(bbox, to_double(line->a()),
        to_double(line->b()), to_double(line->c()));
}

inline bool do_intersect_line_2(
    const Bbox_2 &box, double line_a, double line_b, double line_c);
{
    Bbox_2_Line_2_pair pair(box, line_a, line_b, line_c);
    return pair.intersection_type() != Bbox_2_Line_2_pair::NO;
}

template <class Line>
bool do_intersect_line_2(
    Bbox_2 const &bbox, Line const &line)
{
    return do_intersect_line_2(bbox, to_double(line->a()),
        to_double(line->b()), to_double(line->c()));
}

template <class Line>
bool do_intersect_line_2(
    Line const &line, Bbox_2 const &bbox)
{
    return do_intersect_line_2(bbox, to_double(line->a()),
        to_double(line->b()), to_double(line->c()));
}

template <class R>
inline bool do_intersect(
    const Line_2<R> &line,
    const Bbox_2 &box)
{
    return do_intersect(box, line);
}

CGAL_END_NAMESPACE

#endif
