// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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

    Faces faces;
    Edges boundary_edges;
  };

  Vertex_handle insert_impl(const Point& p, Zone& zone)
  {
#ifdef DEBUG
    std::cerr << "insert(" << p << "): " 
              << zone.boundary_edges.size() << " edges." << std::endl;
#endif
    return triangulation_ref_impl().
      star_hole(p,
		zone.boundary_edges.begin(),
		zone.boundary_edges.end(),
		zone.faces.begin(),
		zone.faces.end());
  }
}; // end Triangulation_mesher_level_traits_2

}; // end namespace CGAL

#endif // CGAL_MESH_2_TRIANGULATION_MESHER_LEVEL_TRAITS_2_H
