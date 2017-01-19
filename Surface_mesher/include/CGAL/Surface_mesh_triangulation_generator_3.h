// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_TRIANGULATION_GENERATOR_3_H
#define CGAL_SURFACE_MESH_TRIANGULATION_GENERATOR_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>

namespace CGAL {

template <class Kernel>
class Surface_mesh_triangulation_generator_3
{
  typedef CGAL::Surface_mesh_vertex_base_3<Kernel> Vb;
  typedef CGAL::Surface_mesh_cell_base_3<Kernel> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
public:
  typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Type;
  typedef Type type; // Boost meta-programming compatibility
};

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_TRIANGULATION_GENERATOR_3_H
