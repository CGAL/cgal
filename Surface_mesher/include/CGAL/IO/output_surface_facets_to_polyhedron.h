// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez

#ifndef CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H
#define CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/value_type_traits.h>

namespace CGAL {

/// \deprecated Gets reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.
///
/// This variant exports the surface as a polyhedron.
/// It requires the surface to be manifold. For this purpose,
/// you may call make_surface_mesh() with Manifold_tag or Manifold_with_boundary_tag parameter.
///
/// @commentheading Template Parameters:
/// @param SurfaceMeshComplex_2InTriangulation_3 model of the SurfaceMeshComplex_2InTriangulation_3 concept.
/// @param Polyhedron an instance of CGAL::Polyhedron_3<Traits>.
///
/// @return true if the surface is manifold and orientable.
template <class SurfaceMeshComplex_2InTriangulation_3,
          class Polyhedron>
CGAL_DEPRECATED_MSG(
        "Please use facets_in_complex_2_to_triangle_mesh() instead."
        ) bool
output_surface_facets_to_polyhedron(
  const SurfaceMeshComplex_2InTriangulation_3& c2t3, ///< Input surface.
  Polyhedron& output_polyhedron) ///< Output polyhedron.
{
  facets_in_complex_2_to_triangle_mesh(c2t3, output_polyhedron);


  // TODO: return true if the surface is manifold and orientable
  return true;
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H
