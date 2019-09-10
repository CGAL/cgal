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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_FILE_TETGEN_H
#define CGAL_IO_FILE_TETGEN_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/IO/File_medit.h>
#include <iostream>
#include <map>
#include <string>
#include <CGAL/utility.h>

namespace CGAL {

namespace Mesh_3 {

template <class C3T3, bool rebind, bool no_patch>
void
output_to_tetgen(std::string filename,
                 const C3T3& c3t3)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "Output to tetgen:\n";
#endif

  typedef Medit_pmap_generator<C3T3,rebind,no_patch> Generator;
  typedef typename Generator::Cell_pmap Cell_pmap;
  typedef typename Generator::Facet_pmap Facet_pmap;
  typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
  typedef typename Generator::Vertex_pmap Vertex_pmap;

  Cell_pmap cell_pmap(c3t3);
  Facet_pmap facet_pmap(c3t3,cell_pmap);
  Facet_pmap_twice facet_pmap_twice(c3t3,cell_pmap);
  Vertex_pmap vertex_pmap(c3t3,cell_pmap,facet_pmap);

  output_to_tetgen(filename,
                   c3t3,
                   vertex_pmap,
                   facet_pmap,
                   cell_pmap,
                   facet_pmap_twice,
                   Generator().print_twice());

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "done.\n";
#endif
}



template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void
output_to_tetgen(std::string filename,
                 const C3T3& c3t3,
                 const Vertex_index_property_map& /* vertex_pmap */,
                 const Facet_index_property_map& /* facet_pmap */,
                 const Cell_index_property_map& /* cell_pmap */,
                 const Facet_index_property_map_twice& /* facet_twice_pmap */ = Facet_index_property_map_twice(),
                 const bool /* print_each_facet_twice */ = false)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Facet Facet;

  const Tr& tr = c3t3.triangulation();

  std::map<Vertex_handle, std::size_t> V;

  //-------------------------------------------------------
  // File output
  //-------------------------------------------------------

  //-------------------------------------------------------
  // nodes
  //-------------------------------------------------------

  std::string node_filename = filename + ".node";
  std::ofstream node_stream(node_filename.c_str());

  node_stream << std::setprecision(17);
  node_stream << tr.number_of_vertices() << " 3 0 0" << std::endl;

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

    node_stream << vert_counter << " " << x << " " << y << " " << z;
    node_stream << std::endl;
  }
  node_stream.close();


  //-------------------------------------------------------
  // Elements
  //-------------------------------------------------------

  std::string elem_filename = filename + ".elem";
  std::ofstream elem_stream(elem_filename.c_str());

  elem_stream << std::setprecision(17);
  elem_stream << c3t3.number_of_cells_in_complex() << " 4 0" << std::endl;

  std::size_t cell_counter = 0;
  for (Cell_iterator
         cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end();
       cit != end; ++cit)
  {
    const Cell_handle ch = cit;

    elem_stream << ++cell_counter;
    for (int i=3; i>=0; i--)
    {
      const Vertex_handle vh = ch->vertex(i);
      elem_stream << " " << V[vh];
    }
    elem_stream << std::endl;
  }
  elem_stream.close();


  //-------------------------------------------------------
  // Face
  //-------------------------------------------------------

  std::string face_filename = filename + ".face";
  std::ofstream face_stream(face_filename.c_str());

  face_stream << std::setprecision(17);
  face_stream << c3t3.number_of_facets_in_complex() << " 0" << std::endl;

  std::size_t facet_counter = 0;
  for(Facet_iterator
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
      fit != end; ++fit )
  {
    const Facet& facet = *fit;

    Vertex_handle vh1 = facet.first->vertex((facet.second+1)%4);
    Vertex_handle vh2 = facet.first->vertex((facet.second+2)%4);
    Vertex_handle vh3 = facet.first->vertex((facet.second+3)%4);

    face_stream << ++facet_counter << " " << V[vh1] << " " << V[vh2] << " " << V[vh3] << std::endl;
  }
  face_stream.close();

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
} // end output_to_tetgen(...)

} // end namespace Mesh_3




/**
 * @brief outputs mesh to tetgen format
 * @param os the stream
 * @param c3t3 the mesh
 * @param rebind if true, labels of cells are rebinded into [1..nb_of_labels]
 * @param show_patches if true, patches are labeled with different labels than
 * cells. If false, each surface facet is written twice, using label of
 * each adjacent cell.
 */
template <class C3T3>
void
output_to_tetgen(std::string filename,
                 const C3T3& c3t3,
                 bool rebind = false,
                 bool show_patches = false)
{
  if ( rebind )
  {
    if ( show_patches )
      Mesh_3::output_to_tetgen<C3T3,true,false>(filename,c3t3);
    else
      Mesh_3::output_to_tetgen<C3T3,true,true>(filename,c3t3);
  }
  else
  {
    if ( show_patches )
      Mesh_3::output_to_tetgen<C3T3,false,false>(filename,c3t3);
    else
      Mesh_3::output_to_tetgen<C3T3,false,true>(filename,c3t3);
  }
}

} // end namespace CGAL

#endif // CGAL_IO_FILE_TETGEN_H
