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

#ifndef CGAL_ZERO_GRADIENT_H
#define CGAL_ZERO_GRADIENT_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup PkgIsosurfacing3Ref
 *
 * \brief Class template for a gradient that is always zero.
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
    /**
     * \ingroup PkgIsosurfacing3Ref
     *
     * \brief Evaluate the gradient at a point in space.
     *
     * \param point the point at which the gradient is computed
     */
    Vector operator()(const Point& point) const {
        return Vector(0, 0, 0);
    }
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_ZERO_GRADIENT_H
