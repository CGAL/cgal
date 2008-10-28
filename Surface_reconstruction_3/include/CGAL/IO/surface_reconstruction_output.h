// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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

#include <CGAL/value_type_traits.h>

CGAL_BEGIN_NAMESPACE


/// Get reconstructed surface out of a SurfaceMeshC2T3 object.
template <class SurfaceMeshC2T3,   ///< Model of SurfaceMeshC2T3.
          typename OutputIterator> ///< OutputIterator value_type must be Triangle_3.
void
output_surface_facets(const SurfaceMeshC2T3& c2t3,
                      OutputIterator output)
{
  typedef typename SurfaceMeshC2T3::Triangulation Tr;
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

      *output = Triangle(points[0],points[1],points[2]);
      output++;
    }
  }
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H
