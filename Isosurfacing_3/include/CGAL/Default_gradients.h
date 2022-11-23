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

#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian_grid_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Template class for a gradient that is always zero.
 *
 * \details This gradient function can be used for Marching Cubes, which does not require a gradient.
 *
 * \tparam GeomTraits the traits for this gradient.
 */
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

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Template class for a gradient that is calculated using finite differences.
 *
 * \details This gradient function evaluates an implicit value function at six points
 *          that are a given distance delta away from the queried point.
 *
 * \tparam GeomTraits the traits for this gradient.
 *
 * \tparam PointFunction the type of the value function
 */
template <class GeomTraits, typename PointFunction>
class Finite_difference_gradient {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

public:
    Finite_difference_gradient(const PointFunction& func, const FT delta = 0.001) : func(func), delta(delta) {}

    Vector operator()(const Point& point) const {  // TODO
        // compute the gradient by sampling the function with finite differences

        const Point p0 = point + Vector(delta, 0, 0);
        const Point p1 = point - Vector(delta, 0, 0);
        const Point p2 = point + Vector(0, delta, 0);
        const Point p3 = point - Vector(0, delta, 0);
        const Point p4 = point + Vector(0, 0, delta);
        const Point p5 = point - Vector(0, 0, delta);

        const FT gx = (func(p0) - func(p1)) / (2 * delta);
        const FT gy = (func(p2) - func(p3)) / (2 * delta);
        const FT gz = (func(p4) - func(p5)) / (2 * delta);

        return Vector(gx, gy, gz);
    }

private:
    const PointFunction func;
    FT delta;
};

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Template class for a gradient that is stored in a Cartesian_grid_3.
 *
 * \details The gradient at any point is calculated using trilinear interpolation.
 *
 * \tparam GeomTraits the traits for this gradient.
 */
template <class GeomTraits>
class Explicit_cartesian_grid_gradient {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point;
    typedef typename Geom_traits::Vector_3 Vector;

    typedef std::shared_ptr<Cartesian_grid_3<Geom_traits>> Grid;

public:
    Explicit_cartesian_grid_gradient(const Grid& grid) : grid(grid) {}

    Vector operator()(const Point& point) const {
        // trilinear interpolation of stored gradients

        const Bbox_3& bbox = grid->get_bbox();
        const Vector& spacing = grid->get_spacing();

        // calculate min index including border case
        std::size_t min_i = (point.x() - bbox.xmin()) / spacing.x();
        std::size_t min_j = (point.y() - bbox.ymin()) / spacing.y();
        std::size_t min_k = (point.z() - bbox.zmin()) / spacing.z();
        if (min_i == grid->xdim() - 1) {
            min_i--;
        }
        if (min_j == grid->ydim() - 1) {
            min_j--;
        }
        if (min_k == grid->zdim() - 1) {
            min_k--;
        }

        const FT min_x = min_i * spacing.x() + bbox.xmin();
        const FT min_y = min_j * spacing.y() + bbox.ymin();
        const FT min_z = min_k * spacing.z() + bbox.zmin();

        const FT f_i = (point.x() - min_x) / spacing.x();
        const FT f_j = (point.y() - min_y) / spacing.y();
        const FT f_k = (point.z() - min_z) / spacing.z();

        const Vector g000 = grid->gradient(min_i + 0, min_j + 0, min_k + 0);
        const Vector g001 = grid->gradient(min_i + 0, min_j + 0, min_k + 1);
        const Vector g010 = grid->gradient(min_i + 0, min_j + 1, min_k + 0);
        const Vector g011 = grid->gradient(min_i + 0, min_j + 1, min_k + 1);
        const Vector g100 = grid->gradient(min_i + 1, min_j + 0, min_k + 0);
        const Vector g101 = grid->gradient(min_i + 1, min_j + 0, min_k + 1);
        const Vector g110 = grid->gradient(min_i + 1, min_j + 1, min_k + 0);
        const Vector g111 = grid->gradient(min_i + 1, min_j + 1, min_k + 1);

        const Vector g0 = g000 * (1 - f_i) * (1 - f_j) * (1 - f_k);
        const Vector g1 = g001 * (1 - f_i) * (1 - f_j) * f_k;
        const Vector g2 = g010 * (1 - f_i) * f_j * (1 - f_k);
        const Vector g3 = g011 * (1 - f_i) * f_j * f_k;
        const Vector g4 = g100 * f_i * (1 - f_j) * (1 - f_k);
        const Vector g5 = g101 * f_i * (1 - f_j) * f_k;
        const Vector g6 = g110 * f_i * f_j * (1 - f_k);
        const Vector g7 = g111 * f_i * f_j * f_k;

        return g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7;
    }

private:
    const Grid grid;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DEFAULT_GRADIENT_H
