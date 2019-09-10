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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SUBDIVIDE_SKIN_SURFACE_MESH_3_H
#define CGAL_SUBDIVIDE_SKIN_SURFACE_MESH_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

// This code is based on the Polyhedron tutorial

template <class SkinSurface_3,
          class Polyhedron_3,
          class SubdivisionPolicy_3>
class Skin_surface_sqrt3
{
  typedef Polyhedron_3                                      Polyhedron;
  typedef SkinSurface_3                                     Skin_surface_3;
  // Projects points to the skin surface:
  typedef SubdivisionPolicy_3                               Subdivision_policy;

  typedef typename Polyhedron::Traits                       Kernel;
  typedef typename Kernel::Point_3                          Point;
  typedef typename Kernel::Vector_3                         Vector;

  typedef typename Polyhedron::Vertex                       Vertex;
  typedef typename Polyhedron::Vertex_handle                Vertex_handle;
  typedef typename Polyhedron::Vertex_iterator              Vertex_iterator;
  typedef typename Polyhedron::Edge_iterator                Edge_iterator;
  typedef typename Polyhedron::Halfedge_handle              Halfedge_handle;
  typedef typename Polyhedron::Halfedge_iterator            Halfedge_iterator;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
                                                            HV_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                                            HF_circulator;
  typedef typename Polyhedron::Facet                        Facet;
  typedef typename Polyhedron::Facet_handle                 Facet_handle;
  typedef typename Polyhedron::Facet_iterator               Facet_iterator;
  typedef typename Kernel::FT                               FT;

public:
  Skin_surface_sqrt3(const SkinSurface_3 &skin,
                     Polyhedron &P,
                     const SubdivisionPolicy_3 &policy)
    : P(P), ss(skin), policy(policy) { }

  //*********************************************
  // Subdivision
  //*********************************************
  int subdivide(int iter=1)
  {
    // check for valid polygon mesh
    if(P.size_of_facets() == 0)
      return false;

    // normalize border: there is no border
    P.normalize_border();

    while (iter > 0) {
      do_subdivide();
      iter--;
    }
    return true;
  }

private:

  //*********************************************
  // Subdivide
  //*********************************************
  void do_subdivide()
  {
    // We use that new vertices/halfedges/facets are appended at the end.
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges
    Facet_iterator last_f = P.facets_end();
    -- last_f;  // the last of the old facets

    // split edges
    Edge_iterator  e = P.edges_begin ();
    do {
      split_halfedge(e);
    } while ( e++ != last_e);

    Vertex_iterator v = P.vertices_begin();
    do {
      Halfedge_handle h_cir, h_start;
      h_cir = h_start = v->halfedge () ;
      do {
        P.split_facet (h_cir->prev(), h_cir->next());
        h_cir = h_cir->next()->opposite();
      } while (h_cir != h_start);
    } while (v++ != last_v);

    v = ++last_v; // First new vertex
    last_v = P.vertices_end();
    do {
      v->point() = policy.to_surface(v);
    } while (++v != last_v);
  }

  //*********************************************
  // Split halfedge
  //*********************************************
  void split_halfedge(Halfedge_handle e)
  {
    // Create a new vertices on e.
    Point p_new = e->vertex()->point();
    e = e->prev();
    p_new = p_new + .5*(e->vertex()->point()-p_new);

    P.split_vertex( e, e->next()->opposite());
    e->next()->vertex()->point() = p_new;
  }

  Polyhedron &P;
  const SkinSurface_3 &ss;
  const Subdivision_policy &policy;
};

template <class SkinSurface_3, class Polyhedron_3>
void subdivide_skin_surface_mesh_3(const SkinSurface_3 &skin,
                                   Polyhedron_3 &p,
                                   int nSubdiv = 1) {
  while (nSubdiv > 0) {
    skin.subdivide_mesh_3(p);
    nSubdiv--;
  }
}

} //namespace CGAL

#endif // CGAL_SUBDIVIDE_SKIN_SURFACE_MESH_3_H
