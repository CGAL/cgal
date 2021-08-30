// Copyright (c) 2003 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_walker.h>

#undef _DEBUG
#define _DEBUG 157
#include <CGAL/Nef_2/debug.h>

typedef CGAL::Gmpz                          NT;
typedef CGAL::Extended_homogeneous_3<NT>    Kernel;
typedef Kernel::Point_3                     Point_3;
typedef Kernel::Segment_3                   Segment_3;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

typedef CGAL::SNC_items<Kernel, bool>       SNC_items;
typedef CGAL::SNC_structure<SNC_items>      SNC_structure;
typedef CGAL::Nef_polyhedron_3<SNC_items>   Nef_polyhedron;
typedef CGAL::SNC_walker<SNC_structure>     SNC_walker;

typedef SNC_structure::Vertex_handle        Vertex_handle;
typedef SNC_structure::Halfedge_handle      Halfedge_handle;
typedef SNC_structure::Halffacet_handle     Halffacet_handle;
typedef SNC_structure::Volume_handle        Volume_handle;
typedef SNC_structure::Object_handle        Object_handle;

int main() {
  CGAL::IO::set_pretty_mode(std::cerr);
  SETDTHREAD(151*157);

  Point_3 source, target;
  std::cin >> source;
  TRACEN( "source:" << source);
  std::cin >> target;
  TRACEN( "target:" << target);

  Point_3 p0( 0, 0, 0, 1);
  Point_3 p1( 4, 0, 0, 1);
  Point_3 p2( 0, 4, 0, 1);
  Point_3 p3( 0, 0, 4, 1);
  Point_3 p4( 4, 4, 4, 1);
  Polyhedron P1, P2;
  P1.make_tetrahedron( p1, p2, p3, p0);
  P2.make_tetrahedron( p3, p2, p1, p4);
  //TRACEN( P1 << P2);

  TRACEN("constructing model...");
  Nef_polyhedron NP1(P1), NP2(P2), NP(NP1^NP2);
  TRACEN("locating source...");
  Object_handle o_src = NP.locate( source);
  TRACEN("walking to target...");
  Object_handle o = NP.locate( o_src, source, target);
  TRACEN("comparing result with ray shooting...");
  CGAL_assertion_code(Object_handle oref = NP.locate(target));

  Vertex_handle v;
  Halfedge_handle e;
  Halffacet_handle f;
  Volume_handle c;
  if( assign( v, o)) { // TODO: enable ostream operator for nef faces
    std::cout << "The point is located in a vertex " <<  std::endl;
    CGAL_assertion_code(Vertex_handle vref);
    CGAL_assertion( assign( vref, oref));
    CGAL_assertion( v == vref);
  }
  else if( assign( e, o)) {
    std::cout << "The point is located in an edge " <<  std::endl;
    CGAL_assertion_code(Halfedge_handle eref);
    CGAL_assertion( assign( eref, oref));
    CGAL_assertion( e == eref || NP.SNCexplorer().twin(eref) == e);
  }
  else if( assign( f, o)) {
    std::cout << "The point is located in a facet " <<  std::endl;
    CGAL_assertion_code(Halffacet_handle fref);
    CGAL_assertion( assign( fref, oref));
    CGAL_assertion( f == fref || NP.SNCexplorer().twin(fref) == f);
  }
  else if( assign( c, o)) {
    std::cout << "The point is located in a volume " << std::endl;
    CGAL_assertion_code( Volume_handle cref);
    CGAL_assertion( assign( cref, oref));
    CGAL_assertion( c == cref);
  }
  else
    CGAL_assertion_msg( 0, "wrong handle");

  return 0;
}


