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
// Author(s)     : Laurent RINEAU


#ifndef CGAL_SURFACE_MESHER_VISITOR_H
#define CGAL_SURFACE_MESHER_VISITOR_H

namespace CGAL {

  namespace Surface_mesher {

    template <
      typename Tr,
      typename Surface_mesher,
      typename Previous_level
      >
    class Visitor {
      Surface_mesher* surface;
      Previous_level* previous;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Previous_level Previous_visitor;

      Visitor(Surface_mesher* surface_,
	      Previous_visitor* p)
        : surface(surface_), previous(p) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      template <class E>
      void before_insertion(E,
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

      Previous_visitor& previous_level()
      {
        return *previous;
      }
      
    }; // end class Visitor

  }  // end namespace Surface_mesher

}  // end namespace CGAL
  
#endif // CGAL_SURFACE_MESHER_VISITOR_H
