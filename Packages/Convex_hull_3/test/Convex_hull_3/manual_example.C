// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

typedef CGAL::Cartesian<double>                                R;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R>  HDS;
typedef CGAL::Polyhedron_default_traits_3<R>                   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits, HDS>                   Polyhedron;
typedef R::Point_3                                       Point_3;

int 
main()
{
  /* generate 250 points randomly on a sphere of radius 100.0 */
  CGAL::Random_points_in_sphere_3< Point_3 > gen(100.0);

  std::vector<Point_3> V;
  CGAL::copy_n( gen, 250, std::back_inserter(V) ); /* copy them to a vector */
  
  Polyhedron P; /* define polyhedron to hold convex hull */

  CGAL::convex_hull_3( V.begin(), V.end(), P);    /* compute convex hull */

  return 0;
}
