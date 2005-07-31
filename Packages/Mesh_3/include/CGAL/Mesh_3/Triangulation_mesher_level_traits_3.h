// Copyright (c) 2004, 2005  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
#define CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H

#include <list>
#include <CGAL/Mesher_level.h>
#include <CGAL/Mesher_level_default_implementations.h>

namespace CGAL {

  namespace Mesh_3 {
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

    } // end namespace Mesh_3::details
  } // end namespace Mesh_3

template <typename Tr>
struct Triangulation_mesher_level_traits_3 :
    public Triangulation_ref_impl<Tr>
{
  typedef Tr Triangulation;

  typedef typename  Mesh_3::details::Type_of_points<typename Tr::Weighted_tag,
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
    typedef std::list<Cell_handle> Cells;
    typedef std::list<Facet> Facets;
  public:
    typedef typename Cells::iterator Cells_iterator;
    typedef typename Facets::iterator Facets_iterator;
    typedef typename Cells::const_iterator Cells_const_iterator;
    typedef typename Facets::const_iterator Facets_const_iterator;

    typedef typename Tr::Locate_type Locate_type;

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
    return triangulation_ref_impl().insert_in_hole(p,
						   zone.cells.begin(),
						   zone.cells.end(),
						   f.first, f.second);
  }

}; // end Triangulation_mesher_level_traits_3

}; // end namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_MESHER_LEVEL_TRAITS_3_H
