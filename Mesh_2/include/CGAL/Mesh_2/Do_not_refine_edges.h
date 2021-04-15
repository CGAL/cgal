// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_2_DO_NOT_REFINE_EDGES_H
#define CGAL_MESH_2_DO_NOT_REFINE_EDGES_H

#include <CGAL/license/Mesh_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_2/Refine_edges.h>

namespace CGAL {

namespace Mesh_2 {

template <
  class Tr,
  class Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
  class Container =
    typename details::Refine_edges_base_types<Tr>::Default_container
>
class Do_not_refine_edges :
    public Refine_edges_base<Tr, Is_locally_conform, Container>
{
  typedef Refine_edges_base<Tr, Is_locally_conform, Container> Super;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Point Point;
  typedef typename Tr::Geom_traits Geom_traits;

  typedef typename Geom_traits::FT FT;

  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Face_circulator Face_circulator;

  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

  using Super::triangulation_ref_impl;

public:
  Do_not_refine_edges(Tr& tr_)
    : Super(tr_) {}

  /** \name FUNCTIONS NEEDED BY Mesher_level OVERIDDEN BY THIS CLASS. */

  void scan_triangulation_impl()
  {
  }

  /**
   * Test if the edges of the boundary are locally conforming.
   * Push which that are not in the list of edges to be conformed.
   */
  Mesher_level_conflict_status
  test_point_conflict_from_superior_impl(const Point& p,
                                         Zone& z)
  {
    if(z.locate_type != Tr::FACE || !z.fh->is_in_domain())
      return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
    for(typename Zone::Edges_iterator eit = z.boundary_edges.begin();
        eit != z.boundary_edges.end(); ++eit)
    {
      const Face_handle& fh = eit->first;
      const int& i = eit->second;

      if(fh->is_constrained(i) && !this->is_locally_conform(this->tr, fh, i, p))
      {
        return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
      }
    }

    return NO_CONFLICT;
  }

}; // end class Do_not_refine_edges

} // end namespace Mesh_2

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_2_DO_NOT_REFINE_EDGES_H
