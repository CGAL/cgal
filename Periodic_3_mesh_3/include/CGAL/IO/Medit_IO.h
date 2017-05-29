// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
#define CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H

namespace CGAL {

namespace Periodic_3_mesh_3 {

namespace IO {

// helper function moving periodic triangles in some canonical positions
// in order to get a surface with fewer "holes"
template<class Triangulation>
typename Triangulation::Periodic_triangle
canonicalize_triangle(const typename Triangulation::Periodic_triangle& pt)
{
  typedef typename Triangulation::Offset            Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  int diffx = std::min(o0.x(), std::min(o1.x(), o2.x()));
  int diffy = std::min(o0.y(), std::min(o1.y(), o2.y()));
  int diffz = std::min(o0.z(), std::min(o1.z(), o2.z()));
  Offset diff_off(diffx, diffy, diffz);

  return CGAL::make_array(std::make_pair(pt[0].first, o0 - diff_off),
                          std::make_pair(pt[1].first, o1 - diff_off),
                          std::make_pair(pt[2].first, o2 - diff_off));
}

template<class Triangulation>
typename Triangulation::Periodic_tetrahedron
canonicalize_tetrahedron(const typename Triangulation::Periodic_tetrahedron& pt)
{
  typedef typename Triangulation::Offset             Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  Offset o3 = pt[3].second;

  int diffx = std::min(std::min(o0.x(), o1.x()), std::min(o2.x(), o3.x()));
  int diffy = std::min(std::min(o0.y(), o1.y()), std::min(o2.y(), o3.y()));
  int diffz = std::min(std::min(o0.z(), o1.z()), std::min(o2.z(), o3.z()));
  Offset diff_off(diffx, diffy, diffz);

  return CGAL::make_array(std::make_pair(pt[0].first, o0 - diff_off),
                          std::make_pair(pt[1].first, o1 - diff_off),
                          std::make_pair(pt[2].first, o2 - diff_off),
                          std::make_pair(pt[3].first, o3 - diff_off));
}

// Writing a restricted Delaunay triangulation to the .mesh file format, which
// can be visualized using medit.
// Writing the triangulation to 8 domains.
template <class Stream, class C3t3>
Stream& write_complex_to_medit(Stream &out, C3t3 &c3t3,
                               unsigned occurence_count = 8)
{
  typedef typename C3t3::Triangulation           Triangulation;
  typedef Triangulation                          Tr;

  typedef typename Triangulation::Iso_cuboid     Iso_cuboid;

  typedef typename Triangulation::Triangle       Triangle;
  typedef typename Triangulation::Tetrahedron    Tetrahedron;

  typedef typename C3t3::Facet_iterator          Facet_iterator;
  typedef typename C3t3::Cell_iterator           Cell_iterator;

  Triangulation& t = c3t3.triangulation();
  int number_of_facets = static_cast<int>(c3t3.number_of_facets());
  int number_of_cells = static_cast<int>(c3t3.number_of_cells());
  int number_of_vertices = 3 * number_of_facets + 4 * number_of_cells;
  out << std::setprecision(17);
  out << "MeshVersionFormatted 1\nDimension 3\nVertices"
      << "\n" << number_of_vertices * occurence_count
      << std::endl;

  Iso_cuboid cb = t.domain();

  for(unsigned j = 0; j < occurence_count; j++ ) {
    for (Facet_iterator it =c3t3.facets_begin(); it!=c3t3.facets_end(); it++) {
      Triangle tri = t.triangle(canonicalize_triangle<Tr>(t.periodic_triangle(*it)));
      for(int i = 0; i < 3; i++) {
        out << tri[i].x() + (j&1) << " "
            << tri[i].y() + ((j&2) >> 1) << " "
            << tri[i].z() + ((j&4) >> 2) << " "
            << 32*j + 1 << std::endl;
      }
    }
  }

  for(unsigned j = 0; j < occurence_count; j++ ) {
    for (Cell_iterator it = c3t3.cells_begin(); it !=c3t3.cells_end(); it++) {
      Tetrahedron tet = t.tetrahedron(canonicalize_tetrahedron<Tr>(t.periodic_tetrahedron( it )));
      for(int i = 0; i < 4; i++) {
        out << tet[i].x() + (j&1) << " "
            << tet[i].y() + ((j&2) >> 1) << " "
            << tet[i].z() + ((j&4) >> 2) << " "
            << 32*j + 1 << std::endl;
      }
    }
  }

  int first_vertex = 1;
  out << "Triangles\n"
      << number_of_facets * occurence_count
      << std::endl;
  const int number_of_vertices_on_facets = number_of_facets * 3;
  for(unsigned j = 0; j < occurence_count; j++ ) {

    for( int i = 0; i < number_of_facets; i++) {
      out << i * 3 + j * number_of_vertices_on_facets + first_vertex  << " "
          << i * 3+1 + j * number_of_vertices_on_facets + first_vertex << " "
          << i * 3+2 + j * number_of_vertices_on_facets + first_vertex << " "
          << 128 * j + 1 << std::endl;
    }
  }

  const int shift = number_of_vertices_on_facets * occurence_count;
  out << "Tetrahedra\n"
      << number_of_cells * occurence_count
      << std::endl;
  const int number_of_vertices_on_cells = number_of_cells * 4;
  for(unsigned j = 0; j < occurence_count; j++ ) {

    Cell_iterator it = c3t3.cells_begin();

    for (int i = 0; i < number_of_cells; i++) {
      out << i * 4 +     j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 1 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 2 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << i * 4 + 3 + j * number_of_vertices_on_cells + first_vertex + shift << " "
          << it->subdomain_index() << std::endl;
      //  << 128 * j + 64 + 1 << std::endl;

      it++;
    }
  }

  out << "0\nEnd";
  return out;
}

} // namespace IO

} // namespace Periodic_3_mesh_3

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
