// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : demo/Convex_hull/ch_quickhull_3_demo.C
// package       : $CGAL_Package: Convex_hull_3 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Convex Hulls and Extreme Points
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: 3D convex hull via quickhull algorithm
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

#include <vector>

// NOTE: the choice of double here for a number type may cause problems 
//       for degenerate point sets
typedef CGAL::Cartesian<double>                   R;
typedef CGAL::Convex_hull_traits_3<R>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;

// define point creator 
typedef CGAL::Point_3<R>                          Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;
typedef CGAL::Random_points_in_sphere_3<Point_3, PointCreator> Generator;

void draw_points_and_hull(const std::vector<Point_3>& points, 
                          const CGAL::Object& object)
{
   std::vector<Point_3>::const_iterator p_it;

   CGAL::Geomview_stream geomview;
   geomview << CGAL::RED;
   for (p_it = points.begin(); p_it != points.end(); p_it++)
   {
      geomview << *p_it;
   }

   CGAL::Segment_3<R>     segment;
   CGAL::Triangle_3<R>    triangle;
   Point_3                point;
   Polyhedron_3           polyhedron;

   geomview << CGAL::BLUE;
   if ( assign(point, object) )
      geomview << point;
   else if ( assign(segment, object) )
      geomview << segment;
   else if ( assign(triangle, object) )
      geomview << triangle;
   else if (assign(polyhedron, object));
      geomview << polyhedron;


   std::cout << "Press any key to end the program: ";
   char wait;
   std::cin >> wait;
}



int main(int argc, char* argv[])
{

  if (argc != 2)
  {
      cerr << "Usage: " << argv[0] << " #points " << endl;
      exit(0);
  }

  int num = atoi(argv[1]);
  if (num < 0) 
  {
     cerr << "Usage: " << argv[0] << " #points " << endl;
     cerr << " #points must be >= 0" << endl;
     exit(0);
  }

  std::vector<Point_3> points;
  Generator gen(100.0);

  // generate num points and copy them to a vector 
  CGAL::copy_n( gen, num, std::back_inserter(points) );
  
  // define object to hold convex hull 
  CGAL::Object ch_object;

  // compute convex hull 
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);
  draw_points_and_hull(points, ch_object);
  return 0;
}
