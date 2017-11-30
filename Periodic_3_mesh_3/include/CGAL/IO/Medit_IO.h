// Copyright (c) 2004-2006, 2009, 2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
#define CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H

#include <CGAL/assertions.h>

#include <boost/unordered_map.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>

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
// By default, 7 copies are used, for a total of 8 instances of the domains.
template <class Stream, class C3t3>
Stream& write_complex_to_medit(Stream& out, C3t3& c3t3,
                               unsigned occurence_count = 8)
{
  // if:
  // -"1": only draws a single instance of the domain
  // -"2": draws 2 occurences of the domain, displaying an additional domain on
  //       the (Ox) axis
  // -"4": draws 4 occurences of the domain, displaying an additional domain on
  //       the (Ox) and (Oy) axes
  // -"8": draws 3 occurences of the domain, displaying an additional domain on
  //       the (Ox), (Oy), and (Oz) axes
  CGAL_precondition(occurence_count == 1 || occurence_count == 2 ||
                    occurence_count == 4 || occurence_count == 8);

  typedef typename C3t3::Triangulation           Triangulation;
  typedef Triangulation                          Tr;

  typedef typename Tr::Bare_point                Bare_point;
  typedef typename Tr::Weighted_point            Weighted_point;

  typedef typename C3t3::Vertex_handle           Vertex_handle;

  typedef typename Tr::Vertex_iterator           Vertex_iterator;
  typedef typename C3t3::Facet_iterator          Facet_iterator;
  typedef typename C3t3::Cell_iterator           Cell_iterator;

  typedef typename Tr::Offset                    Offset;

  const Triangulation& tr = c3t3.triangulation();

  // number of reproductions over the different axes
  int Ox_rn = 1 + (((occurence_count - 1) >> 0) & 1);
  int Oy_rn = 1 + (((occurence_count - 1) >> 1) & 1);
  int Oz_rn = 1 + (((occurence_count - 1) >> 2) & 1);

  int number_of_vertices = static_cast<int>(tr.number_of_vertices());
  int number_of_facets = static_cast<int>(c3t3.number_of_facets());
  int number_of_cells = static_cast<int>(c3t3.number_of_cells());

  out << std::setprecision(17);

  out << "MeshVersionFormatted 1\n"
      << "Dimension 3\n"
      << "Vertices\n"
      << (Ox_rn + 1) * (Oy_rn + 1) * (Oz_rn + 1) * number_of_vertices
      << std::endl;

  // Build the set of points that are needed to draw all canonical elements.

  // On each axis, we repeat n+1 times the point, where 'n' is the number of
  // instances of the mesh that will be printed over that axis. This is because
  // a cell 'c' might have point(c,i) that is equal to v with an offset 2

  boost::unordered_map<Vertex_handle, int> V;
  int inum = 1;

  for(int i=0; i<=Oz_rn; i++)
  {
    for(int j=0; j<=Oy_rn; j++)
    {
      for(int k=0; k<=Ox_rn; k++)
      {
        for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
        {
          if(i == 0 && j == 0 && k == 0)
            V[vit] = inum++;

          const Offset off(k, j, i);
          const Weighted_point& p = vit->point();
          const Bare_point bp = tr.construct_point(p, off);

          int id;
          if(i >= 1 || j >= 1 || k >= 1)
            id = 7;
          else
            id = tr.off_to_int(off);

          out << CGAL::to_double(bp.x()) << ' '
              << CGAL::to_double(bp.y()) << ' '
              << CGAL::to_double(bp.z()) << ' '
              << 32 * id + 1
              << '\n';
        }
      }
    }
  }

  out << "Triangles\n"
      << occurence_count * number_of_facets
      << std::endl;
  for(unsigned j=0; j<occurence_count; j++)
  {
    for(Facet_iterator fit = c3t3.facets_begin(); fit != c3t3.facets_end(); ++fit)
    {
      for (int i=0; i<4; i++)
      {
        if (i == fit->second)
          continue;

        Vertex_handle v = (*fit).first->vertex(i);
        const Offset off = tr.int_to_off(j);
        const Offset combined_off = tr.combine_offsets(
                                      off, tr.int_to_off((*fit).first->offset(i)));
        const int vector_offset = combined_off.x() +
                                  combined_off.y() * (Ox_rn + 1) +
                                  combined_off.z() * (Ox_rn + 1) * (Oy_rn + 1);
        out << vector_offset * number_of_vertices + V[v] << " ";
      }
      out << 128 * j + 1  << '\n';
    }
  }

  out << "Tetrahedra\n"
      << occurence_count * number_of_cells
      << std::endl;
  for(unsigned j=0; j<occurence_count; j++)
  {
    for(Cell_iterator cit = c3t3.cells_begin(); cit !=c3t3.cells_end(); cit++)
    {
      for(int i=0; i<4; ++i)
      {
        const Offset off = tr.int_to_off(j);
        const Offset combined_off = tr.combine_offsets(
                                      off, tr.int_to_off(cit->offset(i)));
        const int vector_offset = combined_off.x() +
                                  combined_off.y() * (Ox_rn + 1) +
                                  combined_off.z() * (Ox_rn + 1) * (Oy_rn + 1);

        out << vector_offset * number_of_vertices + V[cit->vertex(i)] << " ";
      }
      out << cit->subdomain_index() << '\n';
    }
  }

  out << "End" << std::endl;
  return out;
}

} // namespace IO

} // namespace Periodic_3_mesh_3

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
