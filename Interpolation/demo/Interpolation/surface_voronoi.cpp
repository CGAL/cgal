// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Julia Floetotto

#include <CGAL/basic.h>

#ifndef CGAL_USE_GEOMVIEW
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on this platform" << std::endl;
  return 0;
}
#else

#ifdef CGAL_STL_SGI_CC
#define STL_HASH_TABLES
#endif

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Voronoi_intersection_2_traits_3.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

#include <CGAL/squared_distance_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Vector_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>


#include <iostream>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Voronoi_intersection_2_traits_3<K> Gt;
typedef CGAL::Regular_triangulation_2<Gt>        Regular_triangulation;

typedef Regular_triangulation::Face_iterator     Face_iterator;

typedef K::FT                                    Coord_type;
typedef K::Point_3                               Point;
typedef K::Vector_3                              Vector;
typedef K::Segment_3                             Segment;
typedef K::Point_2                               Point_2;
//others
typedef std::vector<Point>                       Point_vector;
typedef std::vector<Point_2>                     Point_2_vector;

typedef CGAL::Aff_transformation_3<K>            Transformation;


//////////////////////
// VISU GEOMVIEW
//////////////////////
template<class Point_vector>
void visu_points(CGAL::Geomview_stream & os, const Point_vector & points)
{
  int n = points.size();
  for(int i=0; i<n ; i++)
    os << points[i];
}

//point generation:
// on a sphere:
void generate_sphere_points(const int& n,
                            const double& r,
                            Point_vector& points,
                            //the test point + normal
                            Point &p, Vector &normal){
  CGAL::Random_points_on_sphere_3<Point> g(r);
  std::copy_n( g, n, std::back_inserter(points));
  p = Point(0,0, r);
  normal = Vector(p - CGAL::ORIGIN);
}

// on a cylinder
void generate_cylinder_points(const int& n,
                              const double& r,
                              const double& height,
                              Point_vector& points,
                              //the test point + normal
                              Point &p, Vector &normal){

  Point_2_vector points_2;
  points_2.reserve(n);
  CGAL::Random_points_on_circle_2<Point_2> g(r);
  std::copy_n( g, n , std::back_inserter(points_2));
  CGAL::Random random;

  double h;
  for(int i=0; i< n; i++){
    h = random.get_double(0.0,height);
    points.push_back(Point(points_2[i].x(),points_2[i].y(),h));
    random.get_double(0.0,height);
  }
  p = Point(r,0,0.5*height);
  normal = Vector(p.x(), p.y(),0);
}

// on a cube
void generate_cube_points(const int& n,
                          const double& r,
                          Point_vector& points,
                          //the test point + normal
                          Point &p, Vector &normal){

  Point_2_vector points_2;
  int m=n/6;
  points_2.reserve(m);
  CGAL::points_on_square_grid_2(r,m,std::back_inserter(points_2),
                                CGAL::Creator_uniform_2<Coord_type,Point_2>());
  //take the tangent plane to the sphere:
  for(int i=0; i < m; i++){
    points.push_back(Point(points_2[i].x(), points_2[i].y(), -r));
    points.push_back(Point(points_2[i].x(), points_2[i].y(), r));
    points.push_back(Point(-r, points_2[i].x(), points_2[i].y()));
    points.push_back(Point(r, points_2[i].x(), points_2[i].y()));
    points.push_back(Point(points_2[i].x(), -r, points_2[i].y()));
    points.push_back(Point(points_2[i].x(), r, points_2[i].y()));
  }
  p = Point(0,0,r);
  normal = Vector(0,0,1);
}



//////////////////////

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  int n=150;
  double r = 5.0;

  // Create test point set. Prepare a vector for 100 points.
  int d;
  std::cout<<"Choose type of surface: 0 -- sphere, 1 -- cylinder,"
           << "2 -- cube "<< std::endl;
  std::cin>> d;
  Point_vector points;
  points.reserve(n);
  //create points
  Point p; Vector normal;
  switch(d){
  case 0: generate_sphere_points(n,r,points, p, normal); break;
  case 1: generate_cylinder_points(n,r, 4*r, points, p, normal); break;
  case 2: generate_cube_points(n,r,points, p, normal); break;
  }


  Transformation transform;
  std::cout <<"Choose type of affine transformation: 0 -- identity,"<<
    " 1 -- rotation, 2 -- translation " << std::endl;
  std::cin >> d;
  switch(d){
  case 0:  transform = Transformation(CGAL::IDENTITY);
    break;
  case 1:
    transform = Transformation(Coord_type(1),Coord_type(0),Coord_type(0),
                               Coord_type(0),
                               Coord_type(0),Coord_type(0.9063),
                               Coord_type(-0.42261826),Coord_type(0),
                               Coord_type(0),Coord_type(0.42261826),
                               Coord_type(0.9063),Coord_type(0));
    break;
  case 2:
    transform = Transformation(CGAL::TRANSLATION, Vector(2,3,-1));
    break;
  }
  //apply affine transformation to points
  Point_vector pts;
  std::transform( points.begin(),points.end(), std::back_inserter(pts),
                  transform);
  points=pts;
  p= transform(p);
  normal =transform(normal);

  //define the triangulation:
  Gt traits = Gt(p,normal);
  Regular_triangulation T(traits);

  //insert the points:
  Gt::Construct_weighted_point_2 p2wp = traits.construct_weighted_point_2_object();

  Point_vector::const_iterator pvit = points.begin(), pvend = points.end();
  for(; pvit!=pvend; ++pvit){
    T.insert(p2wp(*pvit++));
  }

  char ch;
  gv << CGAL::violet();
  visu_points(gv,points);

  gv << CGAL::red() << Segment(p, p+ 0.3*normal);
  gv << CGAL::orange() <<p;

  std::cout << "Visualizing the intersection of "
            << "3D Voronoi diagram with tangent plane at "
            << p << "." << std::endl;
  gv << CGAL::blue();
  T.draw_dual(gv);
  Face_iterator fit = T.finite_faces_begin(),
    fend = T.finite_faces_end();
  for(;fit != fend;fit++)
    gv <<CGAL::black()<<T.dual(fit);

  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 1;
}

#endif // CGAL_USE_GEOMVIEW
