// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT


#ifndef CGAL_SURFACE_MESHER_VISITOR_H
#define CGAL_SURFACE_MESHER_VISITOR_H

#include <CGAL/Surface_mesher/Surface_mesher.h>

namespace CGAL {

  namespace Surface_mesher {

    template <typename Tr, typename Surface_mesher>
    class Visitor {
      Surface_mesher* surface;
      Null_mesh_visitor null_visitor;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Facet Facet;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Null_mesh_visitor Previous_visitor;

      Visitor(Surface_mesher* surface_)
        : surface(surface_) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      void before_insertion(const Facet&,
                            const Point& p,
                            Zone& zone) 
      {
	surface->before_insertion_impl(typename Tr::Facet (), p, zone);
      }

      void after_insertion(const Vertex_handle& v)
      {
	surface->restore_restricted_Delaunay(v);
      }

      template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Null_mesh_visitor& previous_level()
      {
        return null_visitor;
      }
      
    }; // end class Visitor

  }  // end namespace Surface_mesher

}  // end namespace CGAL
  
#endif // CGAL_SURFACE_MESHER_VISITOR_H
