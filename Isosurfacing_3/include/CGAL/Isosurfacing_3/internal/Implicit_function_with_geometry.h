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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
#define CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <memory>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// Wrapper for an implicit function that can only be evaluated at a position and not at a vertex.
// Evaluates the geometry to get the vertex position and passes that to the function.
// Works for all VertexDescriptor types.
template <typename GeomTraits,
          typename Geometry_,
          typename PointFunction>
class Implicit_function_with_geometry
{
public:
  using Geom_traits = GeomTraits;
  using FT = typename Geom_traits::FT;

  using Geometry = std::shared_ptr<Geometry_>;
  using Point_function = std::shared_ptr<PointFunction>;

public:
  // Create a function that uses the geometry to evaluate the function at vertex positions.
  Implicit_function_with_geometry(const Geometry& geom,
                                  const Point_function& func)
    : geom(geom),
      func(func)
  { }

  // Get the value of the function at vertex v
  template <typename VertexDescriptor>
  FT operator()(const VertexDescriptor& v) const
  {
    return func->operator()(geom->operator()(v));
  }

private:
  const Geometry geom;
  const Point_function func;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_FUNCTION_WITH_GEOMETRY_H
