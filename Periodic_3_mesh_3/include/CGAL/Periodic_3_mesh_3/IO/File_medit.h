// Copyright (c) 2004-2006, 2009, 2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mikhail Bogdanov
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
#define CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/IO/File_medit.h>

#include <boost/unordered_map.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace CGAL {

namespace Periodic_3_mesh_3 {

namespace internal {

// Helper functions moving periodic elements in some canonical positions
// in order to get a surface with fewer "holes". Unused for now.
template<class Triangulation>
typename Triangulation::Periodic_triangle
canonicalize_triangle(const typename Triangulation::Periodic_triangle& pt)
{
  typedef typename Triangulation::Offset            Offset;

  Offset o0 = pt[0].second;
  Offset o1 = pt[1].second;
  Offset o2 = pt[2].second;
  int diffx = (std::min)(o0.x(), (std::min)(o1.x(), o2.x()));
  int diffy = (std::min)(o0.y(), (std::min)(o1.y(), o2.y()));
  int diffz = (std::min)(o0.z(), (std::min)(o1.z(), o2.z()));
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

  int diffx = (std::min)((std::min)(o0.x(), o1.x()), (std::min)(o2.x(), o3.x()));
  int diffy = (std::min)((std::min)(o0.y(), o1.y()), (std::min)(o2.y(), o3.y()));
  int diffz = (std::min)((std::min)(o0.z(), o1.z()), (std::min)(o2.z(), o3.z()));
  Offset diff_off(diffx, diffy, diffz);

  return CGAL::make_array(std::make_pair(pt[0].first, o0 - diff_off),
                          std::make_pair(pt[1].first, o1 - diff_off),
                          std::make_pair(pt[2].first, o2 - diff_off),
                          std::make_pair(pt[3].first, o3 - diff_off));
}

} // namespace internal

template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void output_to_medit(std::ostream& os,
                     const C3T3& c3t3,
                     const int occurrence_count,
                     const bool distinguish_copies,
                     const Vertex_index_property_map& vertex_pmap,
                     const Facet_index_property_map& facet_pmap,
                     const Cell_index_property_map& cell_pmap,
                     const Facet_index_property_map_twice& = Facet_index_property_map_twice())
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cout << "Output to medit;\n"
            << "\toccurrences = " << occurrence_count << "\n"
            << "\tdistinguish_copies = " << distinguish_copies << std::endl;
#endif

  // if occurrence_count equals:
  // "1" --> only draws a single instance of the domain.
  // "2" --> draws 2 occurrences of the domain, displaying an additional domain
  //         on the (Ox) axis.
  // "4" --> draws 4 occurrences of the domain, displaying an additional domain
  //         on the (Ox) and (Oy) axes.
  // "8" --> draws 3 occurrences of the domain, displaying an additional domain
  //         on the (Ox), (Oy), and (Oz) axes.
  CGAL_precondition(occurrence_count == 1 || occurrence_count == 2 ||
                    occurrence_count == 4 || occurrence_count == 8);

  typedef typename C3T3::Triangulation           Triangulation;
  typedef Triangulation                          Tr;

  typedef typename Tr::Bare_point                Bare_point;
  typedef typename Tr::Weighted_point            Weighted_point;

  typedef typename C3T3::Vertex_handle           Vertex_handle;
  typedef typename C3T3::Facet                   Facet;
  typedef typename C3T3::Cell_handle             Cell_handle;

  typedef typename Tr::Vertex_iterator           Vertex_iterator;

  typedef typename Tr::Offset                    Offset;

  const Triangulation& tr = c3t3.triangulation();

  // number of reproductions over the different axes
  int Ox_rn = 1 + (((occurrence_count - 1) >> 0) & 1);
  int Oy_rn = 1 + (((occurrence_count - 1) >> 1) & 1);
  int Oz_rn = 1 + (((occurrence_count - 1) >> 2) & 1);

  int number_of_vertices = static_cast<int>(tr.number_of_vertices());
  int number_of_facets = static_cast<int>(c3t3.number_of_facets_in_complex());
  int number_of_cells = static_cast<int>(c3t3.number_of_cells_in_complex());

  // Hardcoded values can be passed here to force more copies
  // Ox_rn = 20; Oy_rn = 20; Oz_rn = 1;

  int occ_mult = Ox_rn * Oy_rn * Oz_rn;

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cout << "Outputting mesh to medit... " << std::endl;
  std::cout << "occurrences over each axis: "
            << Ox_rn << " " << Oy_rn << " " << Oz_rn << std::endl;
  std::cout << number_of_vertices << " vertices" << std::endl;
  std::cout << number_of_facets << " facets" << std::endl;
  std::cout << number_of_cells << " cells" << std::endl;
#endif

  os << std::setprecision(17);

  int medit_number_of_vertices = (Ox_rn + 1) * (Oy_rn + 1) * (Oz_rn + 1) * number_of_vertices;

  os << "MeshVersionFormatted 1\n"
      << "Dimension 3\n"
      << "Vertices\n" << medit_number_of_vertices << std::endl;

  // Build the set of points that is needed to draw all the elements.

  // On each axis, we repeat n+1 times the point, where 'n' is the number of
  // instances of the mesh that will be printed over that axis. This is because
  // a cell 'c' might have point(c,i) that lives in the +1 (in x, y, or z) offset

  boost::unordered_map<Vertex_handle, int> V;
  int inum = 1; // '1' because medit ids start at 1

  for(int i=0; i<=Oz_rn; ++i)
  {
    for(int j=0; j<=Oy_rn; ++j)
    {
      for(int k=0; k<=Ox_rn; ++k)
      {
        for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
        {
          if(i == 0 && j == 0 && k == 0)
            V[vit] = inum++;

          const Offset off(k, j, i);
          const Weighted_point& p = tr.point(vit);
          const Bare_point bp = tr.construct_point(p, off);

          int id;
          if(i >= 1 || j >= 1 || k >= 1)
            id = 7;
          else
            id = tr.off_to_int(off);

          os << CGAL::to_double(bp.x()) << ' '
             << CGAL::to_double(bp.y()) << ' '
             << CGAL::to_double(bp.z()) << ' ';

          if(!distinguish_copies || occ_mult == 1)
            os << get(vertex_pmap, vit) << '\n';
          else
            os << 32 * id + 1 << '\n';
        }
      }
    }
  }

  int medit_number_of_triangles = 2 * occ_mult * number_of_facets;

  os << "Triangles\n" << medit_number_of_triangles << std::endl;
  for(auto cit = c3t3.triangulation().cells_begin();
           cit != c3t3.triangulation().cells_end(); ++cit)
  {
    Cell_handle c = cit;
    for(int s=0; s<4; ++s)
    {
      if(!c3t3.is_in_complex(c, s))
        continue;

      // Offsets in x/y/z
      for(int i=0; i<Oz_rn; ++i)
      {
        for(int j=0; j<Oy_rn; ++j)
        {
          for(int k=0; k<Ox_rn; ++k)
          {
            const Offset off(k, j, i);

            for(int vi=1; vi<4; ++vi) // vertices of the facet in complex
            {
              const int pos = (s+vi) % 4;
              const Offset combined_off = tr.combine_offsets(off, tr.int_to_off(c->offset(pos)));
              const int vector_offset = combined_off.x() +
                                        combined_off.y() * (Ox_rn + 1) +
                                        combined_off.z() * (Ox_rn + 1) * (Oy_rn + 1);

              const Vertex_handle v = c->vertex(pos);
              const int id = vector_offset * number_of_vertices + V[v];

              CGAL_assertion(1 <= id && id <= medit_number_of_vertices);
              os << id << " ";
            }

            // For multiple copies, color to distinguish copies rather than to distinguish subdomains
            if(!distinguish_copies || occ_mult == 1)
              os << get(facet_pmap, Facet(c, s)) << '\n';
            else
              os << 1 + k + 3*j + 9*i << '\n';
          }
        }
      }
    }
  }

  os << "Tetrahedra\n" << occ_mult * number_of_cells << std::endl;
  for(int i=0; i<Oz_rn; ++i)
  {
    for(int j=0; j<Oy_rn; ++j)
    {
      for(int k=0; k<Ox_rn; ++k)
      {
        const Offset off(k, j, i);
        for(auto cit = c3t3.cells_in_complex_begin();
                 cit !=c3t3.cells_in_complex_end(); ++cit)
        {
          for(int l=0; l<4; ++l)
          {
            const Offset combined_off = tr.combine_offsets(off, tr.int_to_off(cit->offset(l)));
            const int vector_offset = combined_off.x() +
                                      combined_off.y() * (Ox_rn + 1) +
                                      combined_off.z() * (Ox_rn + 1) * (Oy_rn + 1);

            const int id = vector_offset * number_of_vertices + V[cit->vertex(l)];
            CGAL_assertion(1 <= id && id <= medit_number_of_vertices);
            os << id << " ";
          }

          // For multiple copies, color to distinguish copies rather than to distinguish subdomains
          if(!distinguish_copies || occ_mult == 1)
            os << get(cell_pmap, cit) << '\n';
          else
            os << 1 + k + 3*j + 9*i << '\n';
        }
      }
    }
  }

  os << "End" << std::endl;
}

template <class C3T3, bool rebind>
void output_to_medit(std::ostream& os,
                     const C3T3& c3t3,
                     const int occurrence_count,
                     const bool distinguish_copies)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cout << "Output to medit:\n";
#endif

  CGAL_precondition(c3t3.triangulation().is_1_cover());

  // periodic meshes always print facets twice because the facet in complex
  // might be on the boundary of the domain
  typedef CGAL::SMDS_3::Medit_pmap_generator<C3T3, rebind, true>      Generator;
  typedef typename Generator::Cell_pmap                               Cell_pmap;
  typedef typename Generator::Facet_pmap                              Facet_pmap;
  typedef typename Generator::Facet_pmap_twice                        Facet_pmap_twice;
  typedef typename Generator::Vertex_pmap                             Vertex_pmap;

  Cell_pmap cell_pmap(c3t3);
  Facet_pmap facet_pmap(c3t3, cell_pmap);
  Facet_pmap_twice facet_pmap_twice(c3t3, cell_pmap);
  Vertex_pmap vertex_pmap(c3t3, cell_pmap, facet_pmap);

  Periodic_3_mesh_3::output_to_medit(os, c3t3, occurrence_count, distinguish_copies,
                                     vertex_pmap, facet_pmap, cell_pmap, facet_pmap_twice);

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cout << "done.\n";
#endif
}

} // namespace Periodic_3_mesh_3

namespace IO {

/**
 * \brief outputs a periodic mesh to the .mesh file format, which can be visualized
 *        using medit. By default, 7 copies are used, for a total of 8 instances of the domains.
 * \param os the stream
 * \param c3t3 the mesh
 * \param occurrence_count the number of copies that are printed
 * \param distinguish_copies if set to `true`, each copy is assigned a unique color.
 *                           Otherwise, all domains are drawn with subdomain index-based colors.
 * \param rebind if set to `true`, labels of cells are rebinded into [1..nb_of_labels]
 */
template <class C3T3>
void output_periodic_mesh_to_medit(std::ostream& os,
                                   const C3T3& c3t3,
                                   const int occurrence_count = 8,
                                   const bool distinguish_copies = true,
                                   const bool rebind = false,
                                   const bool /*show_patches*/ = false) // for backward compatibility
{
  if(rebind)
    Periodic_3_mesh_3::output_to_medit<C3T3, true>(os, c3t3, occurrence_count, distinguish_copies);
  else
    Periodic_3_mesh_3::output_to_medit<C3T3, false>(os, c3t3, occurrence_count, distinguish_copies);
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE
using IO::output_periodic_mesh_to_medit;
#endif

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_IO_FILE_MEDIT_H
