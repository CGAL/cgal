// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
#define CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// Wrapper for an implicit function that can only be evaluated at a position and not at a vertex.
// Evaluates the geometry to get the vertex position and passes that to the function.
template <typename Geometry,
          typename PointFunction>
class Implicit_function_with_geometry
{
  using Point_function = PointFunction;

public:
  // creates a function that uses the geometry to evaluate the function at vertex positions.
  Implicit_function_with_geometry(const Geometry& geom,
                                  const Point_function& func)
    : m_geom(geom),
      m_func(func)
  { }

  // gets the value of the function at vertex `v`
  template <typename VertexDescriptor>
  decltype(auto) /*FT*/ operator()(const VertexDescriptor& v) const
  {
    return m_func.operator()(m_geom.operator()(v));
  }

private:
  const Geometry m_geom;
  const Point_function m_func;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
