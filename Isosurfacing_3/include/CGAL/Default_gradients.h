// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_DEFAULT_GRADIENT_H
#define CGAL_DEFAULT_GRADIENT_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Zero_gradient {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

public:
    Vector operator()(const Point& point) const {
        return Vector(0, 0, 0);
    }
};

template <class GeomTraits, typename Function>
class Finite_difference_gradient {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

public:
    Finite_difference_gradient(const Function& func, const FT delta = 0.001) : func(&func), delta(delta) {}

    Vector operator()(const Point& point) const {  // TODO
        const Point p0 = point + Vector(delta, 0, 0);
        const Point p1 = point - Vector(delta, 0, 0);
        const Point p2 = point + Vector(0, delta, 0);
        const Point p3 = point - Vector(0, delta, 0);
        const Point p4 = point + Vector(0, 0, delta);
        const Point p5 = point - Vector(0, 0, delta);

        const FT gx = (func->operator()(p0) - func->operator()(p1)) / (2 * delta);
        const FT gy = (func->operator()(p2) - func->operator()(p3)) / (2 * delta);
        const FT gz = (func->operator()(p4) - func->operator()(p5)) / (2 * delta);

        return Vector(gx, gy, gz);
    }

private:
    const Function* func;
    FT delta;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DEFAULT_GRADIENT_H
