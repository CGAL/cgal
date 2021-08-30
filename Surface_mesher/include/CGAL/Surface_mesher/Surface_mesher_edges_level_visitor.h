// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU


#ifndef CGAL_SURFACE_MESHER_EDGES_LEVEL_VISITOR_H
#define CGAL_SURFACE_MESHER_EDGES_LEVEL_VISITOR_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>

namespace CGAL {

  namespace Surface_mesher {

    template <
      typename Tr,
      typename Surface_mesher,
      typename Previous_level
      >
    class Edges_level_visitor {
      Surface_mesher* surface_mesher;
      Previous_level* previous;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Previous_level Previous_visitor;

      Edges_level_visitor(Surface_mesher* surface_mesher_,
              Previous_visitor* p)
        : surface_mesher(surface_mesher_), previous(p) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      template <class E>
      void before_insertion(E,
                            const Point& p,
                            Zone zone)
      {
        surface_mesher->remove_edges(p, zone);
      }

      void after_insertion(const Vertex_handle& v)
      {
        surface_mesher->after_insertion_impl(v);
      }

      template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Previous_visitor& previous_level()
      {
        return *previous;
      }

    }; // end class Edges_level_visitor

  }  // end namespace Surface_mesher

}  // end namespace CGAL

#endif // CGAL_SURFACE_MESHER_EDGES_LEVEL_VISITOR_H
