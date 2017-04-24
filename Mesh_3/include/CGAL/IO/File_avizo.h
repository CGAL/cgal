// Copyright (c) 2012  GeometryFactory Sarl (France).
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_FILE_AVIZO_H
#define CGAL_IO_FILE_AVIZO_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/IO/File_medit.h>
#include <iostream>
#include <string>
#include <CGAL/utility.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

namespace Mesh_3 {


template <class C3T3>
void
output_to_avizo(std::ostream& os,
                const C3T3& c3t3)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Weighted_point Weighted_point;

  const Tr& tr = c3t3.triangulation();

  CGAL::Unique_hash_map<Vertex_handle, std::size_t> V;

  //-------------------------------------------------------
  // nodes
  //-------------------------------------------------------

  os << std::setprecision(17);
  os << " # Avizo 3D ASCII 2.0\n\n";
  os << "nNodes " << tr.number_of_vertices() << std::endl;
  os << "nTetrahedra " << c3t3.number_of_cells_in_complex() << std::endl;
  os <<  "Parameters {\n"
    "    Materials {\n"
    "        Material3 {\n"
    "            Id 1,\n"
    "            Color 0 0.835294 0.164706\n"
    "        }\n"
    "        Material4 {\n"
    "            Id 2,\n"
    "            Color 0.862745 0.0901961 0.0901961\n"
    "        }\n"
    "        Material5 {\n"
    "            Id 3,\n"
    "            Color 0.94902 0.847059 0.0901961\n"
    "        }\n"
    "        Material6 {\n"
    "            Id 4,\n"
    "            Color 0.8 0.16 0.698646\n"
    "        }\n"
    "        Material7 {\n"
    "            Id 5,\n"
    "            Color 0.494118 0.494118 1\n"
    "        }\n"
    "        Material8 {\n"
    "            Id 6,\n"
    "            Color 0.227451 0.227451 0.968627\n"
    "        }\n"
    "        Material9 {\n"
    "            Id 7,\n"
    "            Color 0.666667 0.666667 0.666667\n"
    "        }\n"
    "    }\n"
    "}\n"
    "Nodes { float[3] Coordinates } @1\n"
    "Tetrahedra { int[4] Nodes } @2\n"
    "TetrahedronData { byte Materials } @3\n"
    "\n"
    "# Data section follows\n"
    "@1\n";

  std::size_t vert_counter = 0;
  for(Finite_vertices_iterator
        vit = tr.finite_vertices_begin(),
        end = tr.finite_vertices_end();
      vit != end; ++vit)
  {
    const Weighted_point& p = vit->point();
    const double x = CGAL::to_double(p.x());
    const double y = CGAL::to_double(p.y());
    const double z = CGAL::to_double(p.z());

    V[vit] = ++vert_counter;

    os << x << " " << y << " " << z << "\n";
  }



  //-------------------------------------------------------
  // Elements
  //-------------------------------------------------------

  os << "\n@2\n";
  for (Cell_iterator
         cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end();
       cit != end; ++cit)
  {
    const Cell_handle ch = cit;
    os << V[ch->vertex(0)] ;
    os << " " << V[ch->vertex(1)] ;
    os << " " << V[ch->vertex(2)] ;
    os << " " << V[ch->vertex(3)]  << "\n";
  }

  os << "\n@3\n";
  for (Cell_iterator
         cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end();
       cit != end; ++cit)
  {
    os << cit->subdomain_index() << "\n";
  }

} // end output_to_avizo(...)

} // end namespace Mesh_3




/**
 * @brief outputs mesh to avizo format
 * @param os the stream
 * @param c3t3 the mesh
 */
template <class C3T3>
void
output_to_avizo(std::ostream& os,
                 const C3T3& c3t3)
{
  Mesh_3::output_to_avizo(os,c3t3);
}

} // end namespace CGAL

#endif // CGAL_IO_FILE_AVIZO_H
