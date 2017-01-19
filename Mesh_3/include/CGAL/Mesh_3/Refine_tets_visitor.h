// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_REFINE_TETS_VISITOR_H
#define CGAL_MESH_3_REFINE_TETS_VISITOR_H

#include <CGAL/license/Mesh_3.h>


namespace CGAL {

namespace Mesh_3 {

  namespace tets {

    template <
      typename Tr,
      typename Refine_tets,
      typename Previous_level>
    class Refine_facets_visitor {
      Refine_tets* refine_tets;
      Previous_level* previous;
      
    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Facet Facet;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Previous_level Previous_visitor;

      Refine_facets_visitor(Refine_tets* refine_tets_,
			    Previous_visitor* p)
        : refine_tets(refine_tets_), previous(p), active_(false)
      {
      }

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      template <typename E, typename P>
      void before_insertion(const E&,
                            const P&,
                            Zone& zone) 
      {
        if ( active_ )
          refine_tets->before_insertion_handle_cells_in_conflict_zone(zone);
      }

      void after_insertion(const Vertex_handle& v)
      {
        if ( active_ )
          refine_tets->update_star(v);
      }

      template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Previous_visitor& previous_level()
      {
        return *previous;
      }
      
      void activate()
      {
        active_=true;
      }
      
      bool is_active() const
      {
        return active_;
      }
      
    private:
      bool active_;
      
    }; // end class Refine_facets_visitor

  } // end namespace Mesh_3::tets

} // end namespace Mesh_3

} // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_TETS_VISITOR_H
