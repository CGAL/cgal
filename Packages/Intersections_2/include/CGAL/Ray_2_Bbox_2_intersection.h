
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

bool do_intersect_ray_2(
    const Bbox_2 &box, double x, double y, double dx, double dy);

template <class Ray>
bool do_intersect_ray_2(
    const Bbox_2 &box,
    const Ray &ray)
{
    double startx = to_double(ray->start().x());
    double starty = to_double(ray->start().y());
    double dx = to_double(ray->direction().to_vector().x());
    double dy = to_double(ray->direction().to_vector().y());
    return do_intersect_ray_2(box, startx, starty, dx, dy);
}

template <class Ray>
inline bool do_intersect_ray_2(
    const Ray &ray,
    const Bbox_2 &box)
{
    return do_intersect_ray_2(box, ray);
}
CGAL_END_NAMESPACE



#endif
