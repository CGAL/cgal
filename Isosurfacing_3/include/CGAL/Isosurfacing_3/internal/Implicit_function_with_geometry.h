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

#ifndef CGAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
#define CGAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Geometry, typename Function>
class Implicit_function_with_geometry {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;

public:
    Implicit_function_with_geometry(const Geometry& geom, const Function& func) : geom(&geom), func(&func) {}

    template <typename VertexDescriptor>
    FT operator()(const VertexDescriptor& v) const {
        return func->operator()(geom->operator()(v));
    }

private:
    const Geometry* geom;
    const Function* func;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
