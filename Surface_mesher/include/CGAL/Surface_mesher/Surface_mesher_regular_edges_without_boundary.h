// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_SURFACE_MESHER_REGULAR_EDGES_WITHOUT_BOUNDARY_H
#define CGAL_SURFACE_MESHER_REGULAR_EDGES_WITHOUT_BOUNDARY_H

#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>

namespace CGAL {

  namespace Surface_mesher {

  template <
    class C2T3,
    class Surface,
    class SurfaceMeshTraits,
    class Criteria
  >
  class Surface_mesher_regular_edges_without_boundary_base
    : public Surface_mesher_regular_edges_base<C2T3, Surface,
                                               SurfaceMeshTraits, Criteria>
  {
    public:
      typedef Surface_mesher_regular_edges_base<C2T3, Surface,
        SurfaceMeshTraits, Criteria> SMREB;
      typedef C2T3 C2t3;
      typedef typename C2T3::Triangulation Tr;
      typedef typename Tr::Point Point;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;

    public:
    Surface_mesher_regular_edges_without_boundary_base(C2T3& c2t3,
                                                       Surface& surface,
                                                       SurfaceMeshTraits mesh_traits,
                                                       Criteria& criteria)
      : SMREB(c2t3, surface, mesh_traits, criteria)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Surface_mesher_regular_edges_without_boundary_base\n";
#endif
    }

    // Initialization function
    void scan_triangulation_impl() {
      SMREB::scan_triangulation_impl(false);
    }

    void after_insertion_impl(const Vertex_handle v) {
      SMREB::after_insertion_impl(v, false);
    }

  };  // end Surface_mesher_regular_edges_without_boundary_base

  }  // end namespace Surface_mesher

}  // end namespace CGAL


#endif // CGAL_SURFACE_MESHER_REGULAR_EDGES_H

