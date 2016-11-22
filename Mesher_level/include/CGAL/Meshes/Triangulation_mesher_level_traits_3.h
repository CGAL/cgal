// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
// $Id$ $Date$
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
#define CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H

#include <vector>
#include <CGAL/Mesher_level.h>
#include <CGAL/Mesher_level_default_implementations.h>
#include <CGAL/tags.h>

namespace CGAL {

  namespace Meshes {
    namespace details {

      template <typename Tag, typename Tr>
      struct Type_of_points
      {
        typedef typename Tr::Point Point;
      };

      template <typename Tr>
      struct Type_of_points<Tag_true, Tr>
      {
        typedef typename Tr::Weighted_point Point;
      };

    } // end namespace Meshes::details
  } // end namespace Meshes

template <typename Tr>
struct Triangulation_mesher_level_traits_3 :
    public Triangulation_ref_impl<Tr>
{
  typedef Tr Triangulation;

  typedef typename  Meshes::details::Type_of_points<typename Tr::Weighted_tag,
                                                    Tr>::Point Point;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;

  using Triangulation_ref_impl<Tr>::triangulation_ref_impl;

  Triangulation_mesher_level_traits_3(Tr& t)
    : Triangulation_ref_impl<Tr>(t)
  {
  }

  class Zone {
    typedef std::vector<Cell_handle> Cells;
    typedef std::vector<Facet> Facets;
  public:
    typedef typename Cells::iterator Cells_iterator;
    typedef typename Facets::iterator Facets_iterator;
    typedef typename Cells::const_iterator Cells_const_iterator;
    typedef typename Facets::const_iterator Facets_const_iterator;

    typedef typename Tr::Locate_type Locate_type;

    Zone() {
      cells.reserve(64);
      boundary_facets.reserve(32);
      internal_facets.reserve(64);
    }

    Locate_type locate_type;
    Cell_handle cell;
    int i, j;

    Cells cells;
    Facets boundary_facets;
    Facets internal_facets;
  };

  Vertex_handle insert_impl(const Point& p, Zone& zone)
  {
    if( zone.locate_type == Tr::VERTEX 
	) return zone.cell->vertex(zone.i);

    const Facet& f = *(zone.boundary_facets.begin());

    const Vertex_handle v = 
      triangulation_ref_impl().insert_in_hole(p,
                                              zone.cells.begin(),
                                              zone.cells.end(),
                                              f.first, f.second);
#ifdef CGAL_MESH_2_DEBUG_INSERTION_RADIUS
#define CGAL_MESH_3_DEBUG_INSERTION_RADIUS
#endif
#ifdef CGAL_MESH_3_DEBUG_INSERTION_RADIUS
    {
    std::vector<Vertex_handle> vertices;

    triangulation_ref_impl().incident_vertices(v, std::back_inserter(vertices));
    
    typedef typename Tr::Geom_traits::FT FT;

    FT sq_insertion_radius = std::numeric_limits<FT>::infinity();

    for(typename std::vector<Vertex_handle>::const_iterator vit =
          vertices.begin();
        vit != vertices.end();
        ++vit)
      sq_insertion_radius = (CGAL::min)(sq_insertion_radius, 
					CGAL::squared_distance(v->point(),
							       (*vit)->point()) );
    std::cerr << "insertion radius: " << CGAL::sqrt(sq_insertion_radius);
#ifdef CGAL_MESH_3_DIRTY_DEBUG_SPHERES
      std::cerr << " \t\tdistance: " 
                << CGAL::sqrt(CGAL::squared_distance(v->point(), 
                                      typename Tr::Geom_traits::Point_3(CGAL::ORIGIN)));
#endif
    std::cerr << std::endl;
    }
#endif // CGAL_MESH_3_DEBUG_INSERTION_RADIUS

    return v;
  }

}; // end Triangulation_mesher_level_traits_3

} // end namespace CGAL

#endif // CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
