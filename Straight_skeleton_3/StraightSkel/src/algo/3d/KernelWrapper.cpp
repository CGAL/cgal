// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   algo/3d/KernelWrapper.cpp
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#include "algo/3d/KernelWrapper.h"

namespace algo { namespace _3d {

KernelWrapper::KernelWrapper() {
    // intentionally does nothing
}

KernelWrapper::~KernelWrapper() {
    // intentionally does nothing
}

Point3SPtr KernelWrapper::intersection(Plane3SPtr plane1, Plane3SPtr plane2, Plane3SPtr plane3) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane1, *plane2, *plane3);
    if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
        result = KernelFactory::createPoint3(*ipoint);
    } else {
        // Let's investigate...
        if (const CGAL::Line3 *iline = CGAL::object_cast<CGAL::Line3>(&obj))
            std::cerr << "Intersection of 3 planes is a line" << std::endl;
        else if(const CGAL::Plane3 *iplane = CGAL::object_cast<CGAL::Plane3>(&obj))
            std::cerr << "Intersection of 3 planes is a plane" << std::endl;
        else
            std::cerr << "Intersection of 3 planes is... not?" << std::endl;

        CGAL_warning_msg(false, "intersection of 3 planes failed to produce a point");
    }
#else
    result = Point3SPtr(kernel::intersection(&(*plane1), &(*plane2), &(*plane3)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Line3SPtr KernelWrapper::intersection(Plane3SPtr plane1, Plane3SPtr plane2) {
    Line3SPtr result = Line3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane1, *plane2);
    if (const CGAL::Line3 *iline = CGAL::object_cast<CGAL::Line3>(&obj)) {
        result = KernelFactory::createLine3(*iline);
    } else {
        std::cerr << "Intersection is... not?" << std::endl;
        CGAL_warning_msg(false, "intersection of 2 planes failed to produce a line");
    }
#else
    result = Line3SPtr(kernel::intersection(&(*plane1), &(*plane2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::intersection(Plane3SPtr plane, Line3SPtr line) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane, *line);
    if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
        result = KernelFactory::createPoint3(*ipoint);
    } else {
        // Let's investigate
        if (const CGAL::Line3 *iline = CGAL::object_cast<CGAL::Line3>(&obj))
            std::cerr << "Intersection of plane and line is the line itself" << std::endl;
        else
            std::cerr << "Intersection of plane and line is... not?" << std::endl;

        CGAL_warning_msg(false, "intersection of plane and line failed to produce a point");
    }
#else
    result = Point3SPtr(kernel::intersection(&(*plane), &(*line)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::intersection(Sphere3SPtr sphere, Line3SPtr line) {
    Point3SPtr result = Point3SPtr();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    CGAL::FT radius;
#ifdef USE_CGAL
    radius = CGAL::disallowed_sqrt(sphere->squared_radius());
#else
    radius = sphere->getRadius();
#endif
    Vector3SPtr dir = normalize(KernelFactory::createVector3(line));
    Plane3SPtr plane = KernelFactory::createPlane3(p_center, dir);
    Point3SPtr p_intersect = intersection(plane, line);
    CGAL::FT dist = distance(p_center, p_intersect);
    if (dist == radius) {
        result = p_intersect;
    } else if (dist < radius) {
        CGAL::FT amount = - CGAL::disallowed_sqrt(radius*radius - dist*dist);
        result = KernelFactory::createPoint3((*p_intersect) + ((*dir)*amount));
    }
    //DEBUG_SPTR(result);
    return result;
}

Plane3SPtr KernelWrapper::bisector(Plane3SPtr plane1, Plane3SPtr plane2) {
    Plane3SPtr result = Plane3SPtr();
    std::cout << "WARNING: YOU ARE CALLING A BISECTOR, THAT'S AN APPROXIMATE SQRT" << std::endl;
    // CGAL_assertion(false);
#ifdef USE_CGAL
    result = KernelFactory::createPlane3(CGAL::bisector(*plane1, *plane2));
#else
    result = Plane3SPtr(kernel::bisector(&(*plane1), &(*plane2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

CGAL::FT KernelWrapper::squared_distance(Point3SPtr p1, Point3SPtr p2) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::squared_distance(*p1, *p2);
#else
    result = kernel::squared_distance(&(*p1), &(*p2));
#endif
    return result;
}

CGAL::FT KernelWrapper::squared_distance(Segment3SPtr segment, Point3SPtr point) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::squared_distance(*segment, *point);
#else
    result = kernel::squared_distance(&(*segment), &(*point));
#endif
    return result;
}

CGAL::FT KernelWrapper::squared_distance(Line3SPtr line, Point3SPtr point) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::squared_distance(*line, *point);
#else
    result = kernel::squared_distance(&(*line), &(*point));
#endif
    return result;
}

CGAL::FT KernelWrapper::squared_distance(Plane3SPtr plane, Point3SPtr point) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::squared_distance(*plane, *point);
#else
    result = kernel::squared_distance(&(*plane), &(*point));
#endif
    return result;
}

CGAL::FT KernelWrapper::distance(Point3SPtr p1, Point3SPtr p2) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::disallowed_sqrt(CGAL::squared_distance(*p1, *p2));
#else
    result = kernel::distance(&(*p1), &(*p2));
#endif
    return result;
}

CGAL::FT KernelWrapper::distance(Plane3SPtr plane, Point3SPtr point) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::disallowed_sqrt(CGAL::squared_distance(*plane, *point));
#else
    result = kernel::distance(&(*plane), &(*point));
#endif
    return result;
}

CGAL::FT KernelWrapper::distance(Line3SPtr line, Point3SPtr point) {
    CGAL::FT result = 0.0;
#ifdef USE_CGAL
    result = CGAL::disallowed_sqrt(CGAL::squared_distance(*line, *point));
#else
    result = kernel::distance(&(*line), &(*point));
#endif
    return result;
}

Plane3SPtr KernelWrapper::opposite(Plane3SPtr plane) {
    Plane3SPtr result = KernelFactory::createPlane3(plane->opposite());
    DEBUG_SPTR(result);
    return result;
}

Line3SPtr KernelWrapper::opposite(Line3SPtr line) {
    Line3SPtr result = KernelFactory::createLine3(line->opposite());
    DEBUG_SPTR(result);
    return result;
}

Vector3SPtr KernelWrapper::normalize(Vector3SPtr v) {
    Vector3SPtr result;
#ifdef USE_CGAL
    result = KernelFactory::createVector3(*v / CGAL::disallowed_sqrt(v->squared_length()));
#else
    result = KernelFactory::createVector3(v->normalize());
#endif
    return result;
}

Plane3SPtr KernelWrapper::offsetPlane(Plane3SPtr plane, CGAL::FT offset)
{
    Plane3SPtr result = Plane3SPtr();
// #define CGAL_SS3_OLD_CODE_OFFSET_PLANE
#ifndef CGAL_SS3_OLD_CODE_OFFSET_PLANE // this assumes a plane with normalized coefficients (i.e., its normal is normalized)
# ifdef USE_CGAL
    const CGAL::FT a = plane->a();
    const CGAL::FT b = plane->b();
    const CGAL::FT c = plane->c();
    const CGAL::FT d = plane->d();
# else
    const double a = plane->getA();
    const double b = plane->getB();
    const double c = plane->getC();
    const double d = plane->getD();
# endif
    CGAL_assertion_code(CGAL::FT sq_n = a*a + b*b + c*c;)
    // std::cout << "normalization check (post offset): " << sq_n
    //           << " (" << CGAL::to_interval(sq_n).first << "; "
    //           << CGAL::to_interval(sq_n).second << ")" << std::endl;
    // std::cout << "a: " << a << std::endl;
    // std::cout << "b: " << b << std::endl;
    // std::cout << "c: " << c << std::endl;
    // std::cout << "d: " << d << std::endl;

    // inaccuracies during normalization since the sqrt is (usually) not exact
    CGAL_assertion((sq_n - 1) < 1e-5);

    result = KernelFactory::createPlane3(a, b, c, d - offset);
#else // CGAL_SS3_OLD_CODE_OFFSET_PLANE
    const Point3& p = plane->point();
# ifdef USE_CGAL
    Vector3 v_norm = plane->orthogonal_vector();
    Vector3 v_normal = v_norm / CGAL::disallowed_sqrt(v_norm.squared_length());
# else
    Vector3 v_normal = plane->normal().normalize();
# endif
    Point3 p_trans = p + (v_normal * offset);
    Plane3 plane_trans(p_trans, v_normal);
    result = KernelFactory::createPlane3(plane_trans);
#endif // CGAL_SS3_OLD_CODE_OFFSET_PLANE

    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::intersectionOffsetPlanes(Plane3SPtr plane_0,
                                                   CGAL::FT w0,
                                                   Plane3SPtr plane_1,
                                                   CGAL::FT w1,
                                                   Plane3SPtr plane_2,
                                                   CGAL::FT w2,
                                                   Plane3SPtr plane_3,
                                                   CGAL::FT w3)
{
    Point3SPtr result;

    CGAL::FT a0 = plane_0->a(); CGAL::FT b0 = plane_0->b(); CGAL::FT c0 = plane_0->c(); CGAL::FT d0 = plane_0->d();
    CGAL::FT a1 = plane_1->a(); CGAL::FT b1 = plane_1->b(); CGAL::FT c1 = plane_1->c(); CGAL::FT d1 = plane_1->d();
    CGAL::FT a2 = plane_2->a(); CGAL::FT b2 = plane_2->b(); CGAL::FT c2 = plane_2->c(); CGAL::FT d2 = plane_2->d();
    CGAL::FT a3 = plane_3->a(); CGAL::FT b3 = plane_3->b(); CGAL::FT c3 = plane_3->c(); CGAL::FT d3 = plane_3->d();

// #define CGAL_SS3_DEBUG_PLANES_INTERSECTION
#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                  << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                  << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                  << a3 << " " << b3 << " " << c3 << " " << d3 << std::endl;
    std::cout << "Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3 << std::endl;

    std::cout << "CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                   a1, b1, c1, d1,
                                                   a2, b2, c2, d2,
                                                   a3, b3, c3, d3) << std::endl;
#endif

    CGAL_assertion((a0*a0 + b0*b0 + c0*c0 - 1) <= 1e-5);
    CGAL_assertion((a1*a1 + b1*b1 + c1*c1 - 1) <= 1e-5);
    CGAL_assertion((a2*a2 + b2*b2 + c2*c2 - 1) <= 1e-5);
    CGAL_assertion((a3*a3 + b3*b3 + c3*c3 - 1) <= 1e-5);

    CGAL::FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    // @todo, perturbation ensures that this is correct, but in theory we need
    // to handle the degenerate configuration
    // Disabled because it's a runtime speed up
    //
    // if(CGAL::is_zero(den))
    // {
    //   std::cerr << "Warning: degenerate case (coplanar planes)" << std::endl;
    //   return Point3SPtr();
    // }

    // warning: only valid for normalized coefficients!!
    CGAL::FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    CGAL::FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    CGAL::FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "CHECK x|y|z " << x << " " << y << " " << z << std::endl;
#endif

    result = KernelFactory::createPoint3(x, y, z);
    return result;
}

// @todo avoid the duplication (need to avoid recomputing the denominator)
std::pair<Point3SPtr, CGAL::FT> KernelWrapper::intersectionAndTimeOffsetPlanes(Plane3SPtr plane_0,
                                                                               CGAL::FT w0,
                                                                               Plane3SPtr plane_1,
                                                                               CGAL::FT w1,
                                                                               Plane3SPtr plane_2,
                                                                               CGAL::FT w2,
                                                                               Plane3SPtr plane_3,
                                                                               CGAL::FT w3)
{
    Point3SPtr result;

    CGAL::FT a0 = plane_0->a(); CGAL::FT b0 = plane_0->b(); CGAL::FT c0 = plane_0->c(); CGAL::FT d0 = plane_0->d();
    CGAL::FT a1 = plane_1->a(); CGAL::FT b1 = plane_1->b(); CGAL::FT c1 = plane_1->c(); CGAL::FT d1 = plane_1->d();
    CGAL::FT a2 = plane_2->a(); CGAL::FT b2 = plane_2->b(); CGAL::FT c2 = plane_2->c(); CGAL::FT d2 = plane_2->d();
    CGAL::FT a3 = plane_3->a(); CGAL::FT b3 = plane_3->b(); CGAL::FT c3 = plane_3->c(); CGAL::FT d3 = plane_3->d();

#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                  << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                  << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                  << a3 << " " << b3 << " " << c3 << " " << d3 << std::endl;
    std::cout << "Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3 << std::endl;
    std::cout << "CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                   a1, b1, c1, d1,
                                                   a2, b2, c2, d2,
                                                   a3, b3, c3, d3) << std::endl;
#endif

    CGAL_assertion((a0*a0 + b0*b0 + c0*c0 - 1) <= 1e-5);
    CGAL_assertion((a1*a1 + b1*b1 + c1*c1 - 1) <= 1e-5);
    CGAL_assertion((a2*a2 + b2*b2 + c2*c2 - 1) <= 1e-5);
    CGAL_assertion((a3*a3 + b3*b3 + c3*c3 - 1) <= 1e-5);

    CGAL::FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    // @todo perturbation ensures that this is correct, but in theory we need
    // to handle the degenerate configuration
    // Disabled because it's a runtime speed up
    //
    // if(CGAL::is_zero(den))
    // {
    //   std::cerr << "Warning: degenerate case (coplanar planes)" << std::endl;
    //   return Point3SPtr();
    // }

    // warning: only valid for normalized coefficients!!
    CGAL::FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    CGAL::FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    CGAL::FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

    result = KernelFactory::createPoint3(x, y, z);

    CGAL::FT t = (-a0*b1*c2*d3 + a0*b1*c3*d2 + a0*b2*c1*d3 - a0*b2*c3*d1 - a0*b3*c1*d2 + a0*b3*c2*d1 + a1*b0*c2*d3 - a1*b0*c3*d2 - a1*b2*c0*d3 + a1*b2*c3*d0 + a1*b3*c0*d2 - a1*b3*c2*d0 - a2*b0*c1*d3 + a2*b0*c3*d1 + a2*b1*c0*d3 - a2*b1*c3*d0 - a2*b3*c0*d1 + a2*b3*c1*d0 + a3*b0*c1*d2 - a3*b0*c2*d1 - a3*b1*c0*d2 + a3*b1*c2*d0 + a3*b2*c0*d1 - a3*b2*c1*d0) / den;

// #define CGAL_SS3_DEBUG_PLANES_INTERSECTION
#ifdef CGAL_SS3_DEBUG_PLANES_INTERSECTION
    std::cout << "CHECK x|y|z " << x << " " << y << " " << z << std::endl;
    std::cout << a0*x + b0*y + c0*z + d0 - w0*t << std::endl;
    std::cout << a1*x + b1*y + c1*z + d1 - w1*t << std::endl;
    std::cout << a2*x + b2*y + c2*z + d2 - w2*t << std::endl;
    std::cout << a3*x + b3*y + c3*z + d3 - w3*t << std::endl;

    CGAL_assertion(a0*x + b0*y + c0*z + d0 - w0*t == 0);
    CGAL_assertion(a1*x + b1*y + c1*z + d1 - w1*t == 0);
    CGAL_assertion(a2*x + b2*y + c2*z + d2 - w2*t == 0);
    CGAL_assertion(a3*x + b3*y + c3*z + d3 - w3*t == 0);
#endif

    return { result, t };
}

Point3SPtr KernelWrapper::offsetPoint(Point3SPtr point, Vector3SPtr dir, CGAL::FT offset) {
    Point3SPtr result;
#ifdef USE_CGAL
    Vector3 dir_normalized = *dir / CGAL::disallowed_sqrt(dir->squared_length());
    Point3 p_moved = *point + (dir_normalized * offset);
#else
    Point3 p_moved = *point + (dir->normalize() * offset);
#endif
    result = KernelFactory::createPoint3(p_moved);
    DEBUG_SPTR(result);
    return result;
}

Vector3SPtr KernelWrapper::rotateVector(Vector3SPtr vector, Vector3SPtr axis, CGAL::FT angle) {
    Vector3SPtr result;
    Vector3SPtr v_n = KernelWrapper::normalize(axis);
    CGAL::FT n[3];
    for (unsigned int i = 0; i < 3; i++) {
        n[i] = (*v_n)[i];
    }

    CGAL::FT cos_angle = std::cos(CGAL::to_double(angle));
    CGAL::FT sin_angle = std::sin(CGAL::to_double(angle));

    CGAL::FT rotation[3][3];   // http://de.wikipedia.org/wiki/Drehmatrix
    rotation[0][0] = n[0]*n[0] * (1.0-cos_angle) + cos_angle;
    rotation[0][1] = n[0]*n[1] * (1.0-cos_angle) - n[2] * sin_angle;
    rotation[0][2] = n[0]*n[2] * (1.0-cos_angle) + n[1] * sin_angle;
    rotation[1][0] = n[1]*n[0] * (1.0-cos_angle) + n[2] * sin_angle;
    rotation[1][1] = n[1]*n[1] * (1.0-cos_angle) + cos_angle;
    rotation[1][2] = n[1]*n[2] * (1.0-cos_angle) - n[0] * sin_angle;
    rotation[2][0] = n[2]*n[0] * (1.0-cos_angle) - n[1] * sin_angle;
    rotation[2][1] = n[2]*n[1] * (1.0-cos_angle) + n[0] * sin_angle;
    rotation[2][2] = n[2]*n[2] * (1.0-cos_angle) + cos_angle;

    CGAL::FT rotated[3];
    for (unsigned int r = 0; r < 3; r++) {
        rotated[r] = 0.0;
        for (unsigned int c = 0; c < 3; c++) {
            rotated[r] += rotation[r][c] * (*vector)[c];
        }
    }
    result = KernelFactory::createVector3(rotated[0], rotated[1], rotated[2]);
    DEBUG_SPTR(result);
    return result;
}

Plane3SPtr KernelWrapper::rotatePlane(Plane3SPtr plane, Line3SPtr line, CGAL::FT angle) {
    Plane3SPtr result;
    Point3SPtr point;
#ifdef USE_CGAL
    point = KernelFactory::createPoint3(line->point(0));
#else
    point = KernelFactory::createPoint3(line->point());
#endif
    Vector3SPtr dir = KernelFactory::createVector3(line);
    Vector3SPtr normal = KernelFactory::createVector3(plane);
    Vector3SPtr normal_rotated = rotateVector(normal, dir, angle);
    result = KernelFactory::createPlane3(point, normal_rotated);
    DEBUG_SPTR(result);
    return result;
}

int KernelWrapper::side(Plane3SPtr plane, Point3SPtr point) {
    int result = 0;
#ifdef USE_CGAL
    CGAL::Oriented_side side = plane->oriented_side(*point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
#else
    result = plane->side(*point);
#endif
    return result;
}

int KernelWrapper::orientation(Line3SPtr line1, Line3SPtr line2) {
    int result = 0;
#ifdef USE_CGAL
    Vector3 dir1 = line1->to_vector();
    Vector3 dir2 = line2->to_vector();
#else
    Vector3 dir1 = line1->direction();
    Vector3 dir2 = line2->direction();
#endif
    Point3 p0 = line1->point();
    Point3 p1 = p0 + dir1;
    Point3 p2 = line2->point();
    Plane3 plane(p0, p1, p2);
    Point3 point = p2 + dir2;
#ifdef USE_CGAL
    CGAL::Oriented_side side = plane.oriented_side(point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
#else
    result = plane.side(point);
#endif
    return result;
}

double KernelWrapper::angle(Vector3SPtr v1, Vector3SPtr v2) {
    double result = 0.0;
    CGAL::FT arg = 0.0;
#ifdef USE_CGAL
    arg = ((*v1)*(*v2)) / CGAL::disallowed_sqrt(v1->squared_length() * v2->squared_length());
#else
    arg = ((*v1)*(*v2)) / sqrt(v1->squared_length() * v2->squared_length());
#endif
    // fixes issues with floating point precision
    if (arg <= -1.0) {
        result = M_PI;
    } else if (arg >= 1.0) {
        result = 0.0;
    } else {
        result = acos(CGAL::to_double(arg));
    }
    return result;
}

double KernelWrapper::angle(Line3SPtr line1, Line3SPtr line2) {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(line1);
    Vector3SPtr v2 = KernelFactory::createVector3(line2);
    result = angle(v1, v2);
    return result;
}

double KernelWrapper::angle(Plane3SPtr plane, Line3SPtr line) {
    double result = 0.0;
    Vector3SPtr v_plane = KernelFactory::createVector3(plane);
    Vector3SPtr v_line = KernelFactory::createVector3(line);
    result = angle(v_plane, v_line);
    if (result > M_PI/2.0) {
        result = result - M_PI/2.0;
    } else {
        result = M_PI/2.0 - result;
    }
    return result;
}

double KernelWrapper::angle(Plane3SPtr plane1, Plane3SPtr plane2) {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(plane1);
    Vector3SPtr v2 = KernelFactory::createVector3(plane2);
    result = angle(v1, v2);
    return result;
}

bool KernelWrapper::isInside(Point3SPtr p, Point3SPtr p_box_1, Point3SPtr p_box_2) {
    bool result = true;
    for (unsigned int i = 0; i < 3; i++) {
        if ((*p_box_1)[i] < (*p_box_2)[i]) {
            if (!( (*p_box_1)[i] <= (*p)[i] && (*p)[i] <= (*p_box_2)[i] )) {
                result = false;
            }
        } else {
            if (!( (*p_box_2)[i] <= (*p)[i] && (*p)[i] <= (*p_box_1)[i] )) {
                result = false;
            }
        }
    }
    return result;
}

Vector3SPtr KernelWrapper::cross(Vector3SPtr v1, Vector3SPtr v2) {
    Vector3SPtr result = Vector3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createVector3(CGAL::cross_product(*v1, *v2));
#else
    result = KernelFactory::createVector3(v1->cross(*v2));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::projection(Line3SPtr line, Point3SPtr point) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPoint3(line->projection(*point));
#else
    result = Point3SPtr(kernel::projection(&(*line), &(*point)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::projection(Plane3SPtr plane, Point3SPtr point) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPoint3(plane->projection(*point));
#else
    result = Point3SPtr(kernel::projection(&(*plane), &(*point)));
#endif
    DEBUG_SPTR(result);
    return result;
}

int KernelWrapper::comparePoints(Vector3SPtr v_dir, Point3SPtr p_1, Point3SPtr p_2) {
    int result = 0;
    CGAL::FT value = *v_dir * (*p_2 - *p_1);
    if (value > 0.0) {         // angle < M_PI/2.0
        result = -1;
    } else if (value < 0.0) {  // angle > M_PI/2.0
        result = 1;
    }
    return result;
}

Point3SPtr KernelWrapper::replaceCoord(Point3SPtr point, Point3SPtr replacement,
        unsigned int coord) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    if (coord == 0) {
        result = Point3SPtr(new Point3(replacement->x(), point->y(), point->z()));
    } else if (coord == 1) {
        result = Point3SPtr(new Point3(point->x(), replacement->y(), point->z()));
    } else if (coord == 2) {
        result = Point3SPtr(new Point3(point->x(), point->y(), replacement->z()));
    }
#else
    if (coord == 0) {
        result = Point3SPtr(new Point3(replacement->getX(), point->getY(), point->getZ()));
    } else if (coord == 1) {
        result = Point3SPtr(new Point3(point->getX(), replacement->getY(), point->getZ()));
    } else if (coord == 2) {
        result = Point3SPtr(new Point3(point->getX(), point->getY(), replacement->getZ()));
    }
#endif
    DEBUG_SPTR(result);
    return result;
}

} }
