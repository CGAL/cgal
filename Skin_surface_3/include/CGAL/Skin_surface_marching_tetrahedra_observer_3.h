// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SKIN_SURFACE_MARCHING_TETRAHEDRA_OBSERVER_3_H
#define CGAL_SKIN_SURFACE_MARCHING_TETRAHEDRA_OBSERVER_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>

namespace CGAL {

template <class Vertex_iterator,
          class Cell_iterator,
          class Polyhedron_3>
class Skin_surface_marching_tetrahedra_observer_3
  : public Marching_tetrahedra_observer_default_3
             <Vertex_iterator, Cell_iterator, Polyhedron_3>
{
  typedef Polyhedron_3                             Polyhedron;
  typedef Marching_tetrahedra_observer_default_3
    <Vertex_iterator, Cell_iterator, Polyhedron>   Base;

public:
  Skin_surface_marching_tetrahedra_observer_3() : Base() { }
};

template <class Vertex_iterator,
          class Cell_iterator,
          class P_Traits,
          class SkinSurface_3>
class Skin_surface_marching_tetrahedra_observer_3
  <Vertex_iterator,
   Cell_iterator,
   Polyhedron_3<P_Traits,
                Skin_surface_polyhedral_items_3<SkinSurface_3> > >
    : public Marching_tetrahedra_observer_default_3
               <Vertex_iterator,
                Cell_iterator,
                Polyhedron_3<P_Traits,
                             Skin_surface_polyhedral_items_3<SkinSurface_3> > >
{
public:
  typedef Polyhedron_3<P_Traits,
           Skin_surface_polyhedral_items_3<SkinSurface_3> >     Polyhedron;
  typedef Marching_tetrahedra_observer_default_3
    <Vertex_iterator, Cell_iterator, Polyhedron>                Base;

  typedef Cell_iterator                       T_Cell_iterator;
  typedef typename Polyhedron::Vertex_handle  Polyhedron_vertex_handle;
  typedef typename Polyhedron::Facet_handle   Polyhedron_facet_handle;

  Skin_surface_marching_tetrahedra_observer_3() : Base() { }

  Skin_surface_marching_tetrahedra_observer_3(
    const Skin_surface_marching_tetrahedra_observer_3& traits2) : Base(traits2)
  {}

  void after_facet_insertion(T_Cell_iterator ch, Polyhedron_facet_handle fh) {
    fh->tmc_ch = ch;
  }

};

} //namespace CGAL

#endif // CGAL_SKIN_SURFACE_MARCHING_TETRAHEDRA_OBSERVER_3_H
