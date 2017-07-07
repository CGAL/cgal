// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_UNION_OF_BALLS_3_H
#define CGAL_UNION_OF_BALLS_3_H

#include <CGAL/license/Skin_surface_3.h>


#include <CGAL/Skin_surface_base_3.h>
#include <CGAL/triangulate_power_diagram_3.h>

namespace CGAL {

template <class MixedComplexTraits_3>
class Union_of_balls_3
    : public Skin_surface_base_3<MixedComplexTraits_3>
{
  typedef MixedComplexTraits_3                Gt;
  typedef Union_of_balls_3<Gt>                Self;
  typedef Skin_surface_base_3<Gt>             Base;

public:
  typedef MixedComplexTraits_3                 Geometric_traits;

  typedef typename Base::FT                    FT;
  typedef typename Base::Bare_point            Bare_point;
  typedef typename Base::Weighted_point        Weighted_point;
  typedef typename Base::Vector                Vector;

  typedef typename Base::Regular               Regular;
  typedef typename Base::Vertex_handle         Vertex_handle;
  typedef typename Base::Edge                  Edge;
  typedef typename Base::Facet                 Facet;
  typedef typename Base::Facet_circulator      Facet_circulator;
  typedef typename Base::Cell_handle           Cell_handle;
  typedef typename Base::Simplex               Simplex;
  // pair of a del- and vor-simplex
  typedef typename Base::Anchor_point          Anchor_point;

  typedef typename Base::Quadratic_surface     Quadratic_surface;

  typedef typename Base::Vertex_info           Vertex_info;
  typedef typename Base::Cell_info             Cell_info;

  typedef typename Base::TMC                   TMC;
private:
  typedef typename Base::TMC_Vertex_iterator   TMC_Vertex_iterator;
  typedef typename Base::TMC_Cell_iterator     TMC_Cell_iterator;
  typedef typename Base::TMC_Vertex_handle     TMC_Vertex_handle;
  typedef typename Base::TMC_Cell_handle       TMC_Cell_handle;
  typedef typename Base::TMC_Point             TMC_Point;

  using Base::shrink_factor;
  using Base::geometric_traits;
  using Base::regular;
  using Base::triangulated_mixed_complex;

public:
  template < class WP_iterator >
  Union_of_balls_3(WP_iterator begin, WP_iterator end,
                   Gt gt_ = Gt(),
                   bool _verbose = false);

  template <class Polyhedron_3>
  void mesh_skin_surface_3(Polyhedron_3 &p) const {
    Base::mesh_surface_3(p);
  }
  template <class Polyhedron_3>
  void subdivide_skin_surface_mesh_3(Polyhedron_3 &p) const {
    Base::subdivide_mesh_3(p);
  }
};

template <class MixedComplexTraits_3>
template < class WP_iterator >
Union_of_balls_3<MixedComplexTraits_3>::
Union_of_balls_3(WP_iterator begin, WP_iterator end,
                 Gt gt_,
                 bool _verbose)
  : Base(begin, end, 1, false, gt_, _verbose)
{
  // Construct the Triangulated_mixed_complex:
  Triangulated_mixed_complex_observer_3<TMC, Self> observer(shrink_factor());
  triangulate_power_diagram_3(regular(), triangulated_mixed_complex(), observer, _verbose);

  CGAL_assertion(triangulated_mixed_complex().dimension() == 3);
//   { // NGHK: debug code:
//     CGAL_assertion(triangulated_mixed_complex().is_valid());
//     std::vector<TMC_Vertex_handle> ch_vertices;
//     triangulated_mixed_complex().incident_vertices(triangulated_mixed_complex().infinite_vertex(),
//                           std::back_inserter(ch_vertices));
//     for (typename std::vector<TMC_Vertex_handle>::iterator
//            vit = ch_vertices.begin(); vit != ch_vertices.end(); vit++) {
//       CGAL_assertion(sign(*vit) == POSITIVE);
//     }
//   }
}

} //namespace CGAL

#endif // CGAL_UNION_OF_BALLS_3_H
