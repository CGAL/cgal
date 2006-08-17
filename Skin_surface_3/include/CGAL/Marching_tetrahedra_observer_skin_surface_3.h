// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H
#define CGAL_MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H

#include <CGAL/Marching_tetrahedra_observer_default_3.h>

CGAL_BEGIN_NAMESPACE

template <class Vertex_iterator,
	  class Cell_iterator,
	  class Polyhedron_3>
class Marching_tetrahedra_observer_skin_surface_3
  : public Marching_tetrahedra_observer_default_3
  <Vertex_iterator, Cell_iterator, Polyhedron_3>
{
public:
  typedef Marching_tetrahedra_observer_default_3
  <Vertex_iterator, Cell_iterator, Polyhedron_3> Base;

  typedef Polyhedron_3                        Polyhedron;
  
  typedef Cell_iterator                       T_Cell_iterator;
  typedef typename Polyhedron::Vertex_handle  Polyhedron_vertex_handle; 
  typedef typename Polyhedron::Facet_handle   Polyhedron_facet_handle; 

  Marching_tetrahedra_observer_skin_surface_3() : Base() {
  }

  Marching_tetrahedra_observer_skin_surface_3(
    const  Marching_tetrahedra_observer_skin_surface_3& traits2) : Base(traits2)
  {}

  void after_facet_insertion(
    T_Cell_iterator ch,
    Polyhedron_facet_handle fh) {
    fh->sim = ch->mixed_cell();
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_MARCHING_TETRAHEDRA_OBSERVER_SKIN_SURFACE_3_H
