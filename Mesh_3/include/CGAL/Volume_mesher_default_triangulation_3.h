// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_VOLUME_MESHER_DEFAULT_TRIANGULATION_3_H
#define CGAL_VOLUME_MESHER_DEFAULT_TRIANGULATION_3_H

// traits classes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

// regular
#include <CGAL/Regular_triangulation_3.h>

// IO.h must be included before vertex and cell bases.
#include <CGAL/Mesh_3/IO.h>

// vertex and cell bases
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Volume_mesher_cell_base_3.h>

namespace CGAL {
  namespace Mesh_3 {

    class Volume_mesher_default_triangulation_3_generator {

      // traits class
      struct K2 : public CGAL::Exact_predicates_inexact_constructions_kernel {};
      typedef CGAL::Robust_circumcenter_traits_3<K2>  K;
      typedef CGAL::Regular_triangulation_filtered_traits_3<K> Regular_traits;

      // vertex and cell types
      typedef CGAL::Complex_2_in_triangulation_vertex_base_3<Regular_traits> Vb;

      typedef CGAL::Regular_triangulation_cell_base_3<Regular_traits> Cb1;
      typedef CGAL::Surface_mesh_cell_base_3<Regular_traits, Cb1> Cb2;
      typedef CGAL::Volume_mesher_cell_base_3<Regular_traits, Cb2> Cb;

      // triangulation
      typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
      typedef CGAL::Regular_triangulation_3<Regular_traits, Tds> Tr;

    public:
      // result type
      typedef Tr Type;
      typedef Type type;

    }; // end Volume_mesher_default_triangulation_3_generator

  } // end Mesh_3

  typedef Mesh_3::Volume_mesher_default_triangulation_3_generator::Type 
     Volume_mesher_default_triangulation_3;

} // end namespace CGAL  

#endif // CGAL_VOLUME_MESHER_DEFAULT_TRIANGULATION_3_H
