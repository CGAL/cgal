// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/_test_surface_neighbors_3.C
// revision      : 
// revision_date : 

// author(s)     : Julia Floetotto (Julia.Flototto@sophia.inria.fr

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <iostream>
#include <fstream>

#include <utility>

#include <CGAL/surface_neighbor_coordinates_3.h>
#include <CGAL/Voronoi_intersection_2_traits_3.h>

#include <CGAL/point_generators_3.h> 
#include <CGAL/point_generators_2.h> 
#include <CGAL/copy_n.h>
#include <CGAL/function_objects.h>

#include <CGAL/Origin.h>

/////////////////////////////////////////////////////////////
template < class ForwardIterator > 
bool test_norm(ForwardIterator first, ForwardIterator beyond,
	       typename std::iterator_traits<ForwardIterator>::value_type::second_type
	       norm)
{
   typename
     std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
   for(; first !=beyond; first++)
     sum+= first->second;
   
  return(norm==sum); 
}
//////////////////////////////////

template < class ForwardIterator > 
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
		     typename std::iterator_traits<ForwardIterator>::value_type::second_type
		     norm ,
		     const typename
		     std::iterator_traits<ForwardIterator>::value_type::first_type& p,
		     const typename 
		     std::iterator_traits<ForwardIterator>::value_type::second_type&
		     tolerance )
{
  typedef typename
    std::iterator_traits<ForwardIterator>::value_type::first_type Point;
  Point b(0,0,0);
  ForwardIterator it=first;
  for(; it !=beyond; it++)
    b = b+ (it->second/norm) * (it->first - CGAL::ORIGIN);
  return(CGAL::squared_distance(p,b) <= tolerance);  
}



/////////////////////////////////

template < class Triangul>
void
test_coords(const Triangul& T, const  typename
	    Triangul::Geom_traits::Point_3& p, 
	    const typename Triangul::Geom_traits::Vector_3& n, 
	    const typename  Triangul::Geom_traits::FT&
	    tolerance, const int& version )
{
  typedef CGAL::Voronoi_intersection_2_traits_3<typename Triangul::Geom_traits> I_gt;
  //typedef 
  std::vector< std::pair< typename Triangul::Geom_traits::Point_3, 
    typename Triangul::Geom_traits::FT > >  coords;  
  typedef  std::pair< std::back_insert_iterator< 
    std::vector< std::pair< typename Triangul::Geom_traits::Point_3, 
    typename Triangul::Geom_traits::FT > > >, 
    typename Triangul::Geom_traits::FT >                          Result_type;
  
  typename Triangul::Cell_handle start;
  typename Triangul::Geom_traits::FT  norm;
  //test different function calls
  switch(version){
  case 0:{
    Result_type result
      = CGAL::surface_neighbor_coordinates_3(T, p,n,
					     std::back_inserter
					      (coords));
    norm =  result.second; break;}
  case 1: {
    Result_type result  = 
      CGAL::surface_neighbor_coordinates_3(T, p,std::back_inserter
					   (coords),I_gt(p,n)); 
    norm =  result.second; break;}
  //both versions with locate:
  case 2:{
    start = T.locate(p);
    Result_type result  = CGAL::surface_neighbor_coordinates_3(T, p,n,
						   std::back_inserter
						   (coords), start);
    norm =  result.second; break;}
  case 3: { 
    start = T.locate(p);
    Result_type result  = 
      CGAL::surface_neighbor_coordinates_3(T, p,std::back_inserter
					   (coords), I_gt(p,n),start); 
    
    norm =  result.second;break;}
  //taking all points:
  case 4: {
    Result_type result  
      = CGAL::surface_neighbor_coordinates_3(T.points_begin(),
					     T.points_end(),p,n,
					     std::back_inserter
					     (coords),
					     T.geom_traits()); 
    norm =  result.second; break;}
  case 5: {
    Result_type result  
      = CGAL::surface_neighbor_coordinates_3(T.points_begin(),
					     T.points_end(),p,
					     std::back_inserter
					     (coords),I_gt(p,n)); 
    norm =  result.second;
    break;
  }
  default:
    std::cout << "Switch function calls: Nothing is tested. " <<
      std::endl;
  } 
  
  assert(norm>0);
  assert(test_norm( coords.begin(), coords.end(),norm)); 
  assert(test_barycenter(coords.begin(), coords.end(),norm,p, tolerance));
  coords.clear();
}


/////////////////////////////////////////////////////////////


template <class Triangul>
void
_test_surface_neighbors_3_sphere( const Triangul & )
{
  Triangul T;

  int n=200, m=20; 
  double r = 3;

  typedef typename Triangul::Geom_traits          Gt;
  typedef typename Gt::Point_3                    Point;
  
  //test random points in a square of length r:
  std::vector<Point> points;
  points.reserve(n+m);
  
  // Create n+m-4 points on a sphere of radius 2
  CGAL::Random_points_on_sphere_3<Point> g(r);
  CGAL::copy_n( g, n+m, std::back_inserter(points));
  
  for(int i=0; i<n ; i++)
    T.insert(points[i]);
  
  //test with different calls:
  int k=0;
  for(int i=n;i<n+m;i++){
    test_coords(T, points[i],typename
		Gt::Vector_3(points[i]-CGAL::ORIGIN), 
		typename Gt::FT(0.1),++k % 6);
  }
}


////cube case in the interior of a face
template <class Triangul, class Transformation>
void
_test_surface_neighbors_3_cube(const Triangul &, const Transformation&
			       transform, const int n = 75, typename Triangul::Geom_traits::FT
			       tolerance = typename Triangul::Geom_traits::FT(1e-29),
			       bool grid=true)
{
  Triangul T; 

  int m=10; 
  double r = 3;

  typedef typename Triangul::Geom_traits          Gt;
  typedef typename Gt::FT                         Coord_type;
  typedef typename Gt::Point_3                    Point;
  typedef typename Gt::Point_2                    Point_2;
  typedef typename Gt::Vector_3                   Vector;
  
  //generate random points in a square of length r:
  std::vector<Point_2> points_2_data;
  points_2_data.reserve(n);
    	

  if(grid)
    CGAL::points_on_square_grid_2(r,n,std::back_inserter(points_2_data), 
				  CGAL::Creator_uniform_2<Coord_type,Point_2>()); 
  else{
    CGAL::Random_points_in_square_2<Point_2> g(r);
    CGAL::copy_n( g, n, std::back_inserter(points_2_data));
  }

  //generate random points in a square of length r:
  std::vector<Point_2> points_2_test;
  points_2_test.reserve(m);
  CGAL::Random_points_in_square_2<Point_2> g2(r-1.0);
  CGAL::copy_n( g2, m, std::back_inserter(points_2_test));

 
  for(int i=0; i < n; i++){
    T.insert(transform(Point(points_2_data[i].x(), 
			     points_2_data[i].y(), -r)));
    T.insert(transform(Point(points_2_data[i].x(), 
			     points_2_data[i].y(), r)));   
    T.insert(transform(Point(-r, points_2_data[i].x(), 
			     points_2_data[i].y())));
    T.insert(transform(Point(r, points_2_data[i].x(), points_2_data[i].y())));
    T.insert(transform(Point(points_2_data[i].x(), -r, points_2_data[i].y())));
    T.insert(transform(Point(points_2_data[i].x(), r, points_2_data[i].y())));
  }
  int k=0;
  for(int i=0;i<m;i++){
    //test point on z=r plane:
    test_coords(T,transform(Point(points_2_test[i].x(),
				  points_2_test[i].y(), r)),
		transform(Vector(0,0,1)),tolerance, ++k % 6);
    //test point on x=-r plane:
    test_coords(T,transform(Point(-r, points_2_test[i].x(),
				  points_2_test[i].y())),
		transform(Vector(-1,0,0)),tolerance, ++k % 6 );
    //test point on x=r plane:
    test_coords(T,transform(Point(r, points_2_test[i].x(),
				  points_2_test[i].y())),
		transform(Vector(1,0,0)),tolerance,++k % 6 );
    //test point on y=-r plane:
    test_coords(T,transform(Point(points_2_test[i].x(),
				  -r,points_2_test[i].y())),
		transform(Vector(0,-1,0)),tolerance,++k % 6 );
    //test point on y=r plane:
    test_coords(T,transform(Point(points_2_test[i].x(),
				  r,points_2_test[i].y())),
		transform(Vector(0,1,0)),tolerance,++k % 6);
  }

  //test a sample point: 
  //with Delaunay triangulation filering:
  test_coords(T,transform(Point(points_2_data[n/2].x(),
				points_2_data[n/2].y(), r)),
	      transform(Vector(0,0,1)),Coord_type(0),0);
  //considering all points:
  test_coords(T,transform(Point(points_2_data[n/2].x(),
				 points_2_data[n/2].y(), r)),
	       transform(Vector(0,0,1)),Coord_type(0),4);
  
}
//end of file
