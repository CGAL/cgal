
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
// file          : include/CGAL/Ray_2_Bbox_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_RAY_2_BBOX_2_INTERSECTION_H
#define CGAL_RAY_2_BBOX_2_INTERSECTION_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Ray_2.h>
// #include <CGAL/Segment_2.h>
// #include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

class Bbox_2_Ray_2_pair_impl;

class Bbox_2_Ray_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    ~Bbox_2_Ray_2_pair() ;
    Bbox_2_Ray_2_pair() ;
    Bbox_2_Ray_2_pair(Bbox_2_Ray_2_pair const &o) ;
    Bbox_2_Ray_2_pair(Bbox_2 const &box,
                      double x, double y, double dx, double dy) ;
    Bbox_2_Ray_2_pair& operator=(Bbox_2_Ray_2_pair const &o) ;
    Intersection_results intersection_type() const;
    bool intersection(double &x, double &y) const;
    bool intersection(double &x1, double &y1, double &x2, double &y2) const;
protected:
    Bbox_2_Ray_2_pair_impl *pimpl;
};

/*
template <class R>
inline bool do_intersect(
    const Ray_2<R> &ray,
    const Bbox_2 &box)
{
    typedef Ray_2_Bbox_2_pair<R> pair_t;
    pair_t pair(&ray, &box);
    return pair.intersection_type() != pair_t::NO;
}
*/
CGAL_END_NAMESPACE


CGAL_BEGIN_NAMESPACE
/*
template <class R>
class Bbox_2_Ray_2_pair: public Ray_2_Bbox_2_pair<R> {
public:
    Bbox_2_Ray_2_pair() {}
    Bbox_2_Ray_2_pair(Bbox_2 const *box, Ray_2<R> const *ray)
                :Ray_2_Bbox_2_pair<R> (ray, box){}
};

template <class R>
inline bool do_intersect(
    const Bbox_2 &box,
    const Ray_2<R> &ray)
{
    typedef Bbox_2_Ray_2_pair<R> pair_t;
    pair_t pair(&box, &ray);
    return pair.intersection_type() != pair_t::NO;
}
*/
CGAL_END_NAMESPACE

#endif
