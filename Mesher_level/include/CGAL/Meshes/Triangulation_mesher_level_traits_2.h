// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_2_H
#define CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_2_H

#include <list>
#include <CGAL/Mesher_level.h>

namespace CGAL {

template <typename Tr>
struct Triangulation_mesher_level_traits_2 : 
    public Triangulation_ref_impl<Tr>
{
  typedef Tr Triangulation;
  typedef typename Tr::Point Point;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;

  using Triangulation_ref_impl<Tr>::triangulation_ref_impl;

  Triangulation_mesher_level_traits_2(Tr& tr)
    : Triangulation_ref_impl<Tr>(tr)
  {
  }

  struct Zone {
    typedef std::list<Face_handle> Faces;
    typedef std::list<Edge> Edges;
  public:
    typedef typename Faces::iterator Faces_iterator;
    typedef typename Edges::iterator Edges_iterator;

    typename Tr::Locate_type locate_type;
    Face_handle fh;
    int i;
    Face_handle parent_face; // When the refinement point is the
                             // circumcenter of a face, that member stores
                             // the face.
    Faces faces;
    Edges boundary_edges;
  };

  Vertex_handle insert_impl(const Point& p, Zone& zone)
  {
#ifdef CGAL_MESH_2_DEBUG_INSERTIONS
    std::cerr << "insert(" << p << "): " 
              << zone.boundary_edges.size() << " edges." << std::endl;
#endif
    if( zone.locate_type == Tr::VERTEX )
      return zone.fh->vertex(zone.i);
    return triangulation_ref_impl().
      star_hole(p,
		zone.boundary_edges.begin(),
		zone.boundary_edges.end(),
		zone.faces.begin(),
		zone.faces.end());
  }
}; // end Triangulation_mesher_level_traits_2

} // end namespace CGAL

#endif // CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_2_H
