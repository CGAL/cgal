// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_DEFAULT_TRIANGULATION_3_H
#define CGAL_SURFACE_MESH_DEFAULT_TRIANGULATION_3_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

// traits class
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

// Delaunay
#include <CGAL/Delaunay_triangulation_3.h>

// vertex and cell bases
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

namespace CGAL {
  namespace Surface_mesher {

    template <typename K2 = Exact_predicates_inexact_constructions_kernel>
    class Surface_mesh_default_triangulation_3_generator {

      // traits class
      typedef Robust_circumcenter_traits_3<K2>  K;

      // vertex and cell types
      typedef Surface_mesh_vertex_base_3<K> Vb;
      typedef Surface_mesh_cell_base_3<K> Cb;
      typedef Delaunay_triangulation_cell_base_with_circumcenter_3<K, Cb>
        Cb_with_circumcenter;

      // triangulation
      typedef Triangulation_data_structure_3<Vb, Cb_with_circumcenter> Tds;
      typedef Delaunay_triangulation_3<K, Tds> Tr;

    public:
      // result type
      typedef Tr Type;
      typedef Type type;

    }; // end Surface_mesh_default_triangulation_3_generator

  } // end Surface_mesher

  typedef Surface_mesher::Surface_mesh_default_triangulation_3_generator<>::Type
     Surface_mesh_default_triangulation_3;

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_DEFAULT_TRIANGULATION_3_H
