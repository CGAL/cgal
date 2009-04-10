// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H
#define CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h>
#include <CGAL/value_type_traits.h>

CGAL_BEGIN_NAMESPACE


/// Get reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.
///
/// This variant exports the surface as a triangle soup.
///
/// @commentheading Template Parameters:
/// @param SurfaceMeshComplex_2InTriangulation_3 model of the SurfaceMeshComplex_2InTriangulation_3 concept.
/// @param OutputIterator value_type must be Triangle_3.
///
/// @return true on success.
template <class SurfaceMeshComplex_2InTriangulation_3,
          typename OutputIterator>
void
surface_reconstruction_output_surface_facets(
  const SurfaceMeshComplex_2InTriangulation_3& c2t3, ///< Input surface.
  OutputIterator output_iterator) ///< Output iterator.
{
  typedef SurfaceMeshComplex_2InTriangulation_3 C2t3;
  typedef typename C2t3::Triangulation Tr;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename value_type_traits<OutputIterator>::type Triangle;

  const Tr& tr = c2t3.triangulation();

  for(Finite_facets_iterator fit = tr.finite_facets_begin();
    fit != tr.finite_facets_end();
    ++fit)
  {
    if((*fit).first->is_facet_on_surface((*fit).second) == true)
    {
      std::vector<Point> points(3);
      unsigned int index = 0;
      for(int i=0; i<4; i++)
        if(i != (*fit).second)
          points[index++] = (*fit).first->vertex(i)->point();

      *output_iterator = Triangle(points[0],points[1],points[2]);
      output_iterator++;
    }
  }
}


/// Get reconstructed surface out of a SurfaceMeshComplex_2InTriangulation_3 object.
///
/// This variant exports the surface as a polyhedron. 
/// It requires the surface to be manifold. For this purpose,
/// you may call make_surface_mesh() with Manifold_tag or Manifold_with_boundary_tag parameter.
///
/// @commentheading Template Parameters:
/// @param SurfaceMeshComplex_2InTriangulation_3 model of the SurfaceMeshComplex_2InTriangulation_3 concept.
/// @param PolyhedronTraits_3, PolyhedronItems_3, HalfedgeDS, Alloc see Polyhedron_3 declaration.
///
/// @return true on success.
template <class SurfaceMeshComplex_2InTriangulation_3,
          class PolyhedronTraits_3, class PolyhedronItems_3, class HalfedgeDS, class Alloc>
void
surface_reconstruction_output_surface_facets(
  const SurfaceMeshComplex_2InTriangulation_3& c2t3, ///< Input surface.
  Polyhedron_3<PolyhedronTraits_3, PolyhedronItems_3, HalfedgeDS, Alloc>& output_polyhedron) ///< Output polyhedron.
{
  typedef SurfaceMeshComplex_2InTriangulation_3 C2t3;
  typedef Polyhedron_3<PolyhedronTraits_3, PolyhedronItems_3, HalfedgeDS, Alloc> Polyhedron;

  Complex_2_in_triangulation_3_polyhedron_builder<C2t3, Polyhedron>  builder(c2t3);    
  output_polyhedron.delegate(builder);
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H
