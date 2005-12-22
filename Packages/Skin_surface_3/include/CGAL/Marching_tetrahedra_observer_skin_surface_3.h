// Copyright (c) 2005 RuG (Netherlands)
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H
#define MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H

CGAL_BEGIN_NAMESPACE

template <class Triangulation_3,
	  class Polyhedron_3 >
class Marching_tetrahedra_observer_skin_surface_3 {
public:
  typedef Triangulation_3                       Triangulation;
  typedef typename Triangulation::Cell_handle   Triangulation_cell_handle;
  typedef typename Polyhedron_3::Vertex_handle  Polyhedron_vertex_handle; 
  typedef typename Polyhedron_3::Facet_handle   Polyhedron_facet_handle; 

  Marching_tetrahedra_observer_skin_surface_3() {
  }

  Marching_tetrahedra_observer_skin_surface_3(
    const  Marching_tetrahedra_observer_skin_surface_3& traits2) {
  }

  void after_vertex_insertion(
    Triangulation_cell_handle ch,  int i, int j,
    Polyhedron_vertex_handle vh) {
  }
  
  void after_facet_insertion(
    Triangulation_cell_handle ch,
    Polyhedron_facet_handle fh) {
  }

};

CGAL_END_NAMESPACE

#endif // MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H
