// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Pierre Alliez

#ifndef CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H
#define CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>
#include <CGAL/value_type_traits.h>

namespace CGAL {


/// Gets reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.
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
bool
output_surface_facets_to_polyhedron(
  const SurfaceMeshComplex_2InTriangulation_3& c2t3, ///< Input surface.
  Polyhedron& output_polyhedron) ///< Output polyhedron.
{
  typedef SurfaceMeshComplex_2InTriangulation_3 C2t3;

  Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron>  builder(c2t3);
  output_polyhedron.delegate(builder);

  // TODO: return true if the surface is manifold and orientable
  return true;
}


} //namespace CGAL

#endif // CGAL_OUTPUT_SURFACE_FACETS_TO_POLYHEDRON_H
