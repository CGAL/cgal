// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 21
//
// file          : 
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : chull_traits.lw
// revision      : 2.3  
// revision_date : 01 Feb 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CONVEX_HULL_3_H
#define CGAL_CONVEX_HULL_3_H

#include <CGAL/dd_geo/chull_traits_3.h>
#include <CGAL/dd_geo/chull_support_3.h>

CGAL_BEGIN_NAMESPACE
template <class InputIterator, class Polyhedron>
void
convex_hull_3( InputIterator first, InputIterator beyond, Polyhedron& P)
{
  typedef typename Polyhedron::Traits       PolyTraits;
  typedef typename Polyhedron::Halfedge_data_structure
                                            HDS;
  typedef typename PolyTraits::R            R;
  typedef typename PolyTraits::Point        Point;
  typedef chull_traits_3<R>                 ChullTraits;
#ifndef CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS
  typedef chull< ChullTraits >              ChullType;
#else
  typedef Plane_3<R>                        Plane;
  typedef chull< ChullTraits, Point, Plane> ChullType;
#endif // CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS

  ChullType CH(3);
  for ( ; first != beyond ; ++first)  CH.insert(*first);

#ifndef DDGEO_STL_ITERATORS
  GRAPH<Point,int> G; 
  d3_surface_map( CH, G);
  G.compute_faces();
  Build_polyhedron_from_GRAPH< HDS>  polyhedron_maker( G );
#else
  Build_polyhedron_from_chull< HDS, ChullType>  polyhedron_maker( CH);
#endif // DDGEO_STL_ITERATORS
  Polyhedron P_local;
  P_local.delegate( polyhedron_maker );
  P = P_local;
}
CGAL_END_NAMESPACE

#endif // CGAL_CONVEX_HULL_3_H
