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

#ifndef CGAL_SKIN_SURFACE_3_H
#define CGAL_SKIN_SURFACE_3_H

#include <CGAL/Skin_surface_base_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>
#include <CGAL/FPU.h>

namespace CGAL { 

template <class MixedComplexTraits_3> 
class Skin_surface_3 : public Skin_surface_base_3<MixedComplexTraits_3> {
  typedef MixedComplexTraits_3            Gt;
  typedef Skin_surface_3<Gt>              Self;
  typedef Skin_surface_base_3<Gt>         Base; 
public:
  typedef MixedComplexTraits_3            Geometric_traits;

  typedef typename Gt::Weighted_point     Weighted_point;
  typedef typename Gt::Bare_point         Bare_point;
  typedef typename Gt::FT                 FT;
  // For normal():
  typedef typename Gt::Vector_3           Vector;


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

public:
  using Base::shrink_factor;
  using Base::geometric_traits;
  using Base::regular;
  using Base::triangulated_mixed_complex;
public:
  template < class WP_iterator >
  Skin_surface_3(WP_iterator begin, WP_iterator end, 
                 FT shrink,
                 bool grow_balls = true,
                 Gt gt_ = Gt(),
                 bool _verbose = false
                 );

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
Skin_surface_3<MixedComplexTraits_3>::
Skin_surface_3(WP_iterator begin, WP_iterator end, 
               FT shrink,
               bool grow_balls,
               Gt gt_,
               bool verbose_) 
  : Base(begin, end, shrink, grow_balls, gt_, verbose_) {
    
  // Construct the Triangulated_mixed_complex:
  Triangulated_mixed_complex_observer_3<TMC, Self> observer(shrink_factor());
  triangulate_mixed_complex_3(regular(), shrink_factor(), 
                              triangulated_mixed_complex(),
                              observer, verbose_);

  CGAL_assertion(triangulated_mixed_complex().dimension() == 3);
}

} //namespace CGAL

#endif // CGAL_SKIN_SURFACE_3_H
