// ============================================================================
//
// Copyright (c) 1998-1999 The CGAL Consortium
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
// file          : demo_2/
// revision      : $Revision$
// author(s)     : Julia Floetotto
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================
// Geomview doesn't work on M$ at the moment, so we don't compile this
// file.
//**********************
//demo 2D Interpolation over the plane - using 2D natural neighbors
//**********************


#if defined(__BORLANDC__) || defined(_MSC_VER)
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on Windows, so this demo doesn't work"
            << std::endl;
  return 0;
}
#else

#include <CGAL/basic.h>
#include <fstream>
#include <iostream>
#include <utility>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//#include <CGAL/double.h>
//#include <CGAL/Random.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/interpolation_functions.h>

#include <CGAL/point_generators_2.h> 
#include <CGAL/copy_n.h>

#include <CGAL/IO/Geomview_stream.h>


struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef CGAL::Delaunay_triangulation_2<K>   Dt;
typedef CGAL::Interpolation_traits_2<K>     ITraits;

typedef K::FT                               Coord_type;
typedef K::Point_2                          Point_2;
typedef K::Vector_2                         Vector_2;

typedef K::Point_3                          Point_3;
typedef K::Vector_3                         Vector_3;
typedef K::Segment_3                        Segment_3;
typedef K::Triangle_3                       Triangle_3;

typedef std::map<Point_2, Coord_type, K::Less_xy_2>     Point_value_map ;
typedef std::map<Point_2, Vector_2, K::Less_xy_2 >      Point_vector_map;

typedef std::vector<Point_3>                            Point_vector_3;
typedef std::vector<Point_2>                            Point_vector_2;

////////////////////// 
// VISU GEOMVIEW
////////////////////// 

template<class Point_vector>
void visu_points(CGAL::Geomview_stream & gv, const Point_vector & points)
{
  int n = points.size();
  for(int i=0; i<n ; i++)
    gv << points[i];
}
template < class Point_vector>
void
visu_graph(CGAL::Geomview_stream & gv, const Point_vector& points, int m) 
{

  //generate list of triangles:
  std::vector< Triangle_3 > tr;
  tr.reserve(2*(m-1)*(m-1));
  
  //indices
  for(int i=0; i< m-1;i++)
    for(int j=0; j< m-1; j++){
      tr.push_back(Triangle_3(points[i*m + j],
			      points[j+1 + i*m],
			      points[(j+1) + (i+1)*m]));
      tr.push_back(Triangle_3(points[i*m + j],points[(j+1) + (i+1)*m],
			      points[j + (i+1)*m]));
    } 
  gv.draw_triangles(tr.begin(), tr.end());  
}

//////////////////////////////////////////////////////////////////////////
////POINT GENERATION:
void generate_grid_points(Point_vector_2& points, int m, float h) 
{
  
  //int n = (m+1)*(m+1);
  int n = m*m;
  points.clear();
  points.reserve(n);

  std::cout <<" generate " <<n   << " grid points in square of height " 
	    << h <<std::endl;
  // Create n points from a 16 x 16 grid. Note that the double
  // arithmetic _is_ sufficient to produce exact integer grid points.
  // The distance between neighbors is 34 pixel = 510 / 15.
  CGAL::points_on_square_grid_2((double) h, n,
				std::back_inserter(points),
				CGAL::Creator_uniform_2<double,Point_2>());
  
}

////////////////////// 

int main(int argc,  char* argv[])
{

  //parameters:
  int m = 78;
  double h = 10;
  double g = 2;
  Coord_type w(4);
  Point_vector_2 sample;
  //3D function graph: 2D points + function value in z-direction:
  Point_vector_3 sample_3;

  Dt T;
  Point_value_map values;
  Point_vector_map gradients;

  sample.push_back( Point_2(h,h));
  sample.push_back( Point_2( -h,h));
  sample.push_back( Point_2( h,-h ));
  sample.push_back( Point_2( -h,-h));
  sample.push_back( Point_2( 0.0, h/2));
  sample.push_back( Point_2( 0.0,-h/2));
  sample.push_back( Point_2( h/2, 0.0 ));
  sample.push_back( Point_2( -h/2, 0.0));
  sample.push_back( Point_2(0.0,0.0));

  int s = sample.size();
  for(int j=0; j<s ; j++){
    T.insert(sample[j]);
    
    values.insert(std::make_pair(sample[j], Coord_type(0)));
    gradients.insert(std::make_pair(sample[j], CGAL::NULL_VECTOR));
    sample_3.push_back( Point_3(sample[j].x(),sample[j].y(),0));
  }

  Point_2 p = Point_2(h/3,h/3);
  T.insert(p);
  values.insert(std::make_pair(p, w));
  gradients.insert(std::make_pair(p, Vector_2(-g,-g) ));
   
  p = Point_2(-h/3,-h/3);
  T.insert(p);
  values.insert(std::make_pair(p,w));
  gradients.insert(std::make_pair(p, Vector_2(g,g) ));
   
  p = Point_2(-h/3,h/3);
  T.insert(p);
  values.insert(std::make_pair(p, w));
  gradients.insert(std::make_pair(p, Vector_2(g,-g) ));
   
  p = Point_2(h/3,-h/3);
  T.insert(p);
  values.insert(std::make_pair(p, w));
  gradients.insert(std::make_pair(p, Vector_2(-g,g) ));
   
  sample_3.push_back( Point_3(h/3,h/3,w)); 
  sample_3.push_back( Point_3(-h/3,-h/3,w));
  sample_3.push_back( Point_3(-h/3,h/3,w));
  sample_3.push_back( Point_3(h/3,-h/3,w));
 
  //Interpolation grid:
  Point_vector_2 points;
  generate_grid_points(points, m, 0.999* h); 
  Point_vector_3 points_3;
  
  int method;
  std::cout << "Enter the choice of the interpolation function: "<<std::endl;
  std::cout << " 0 -- linear interpolation" <<std::endl;  
  std::cout << " 1 -- simple quadratic interpolant"  <<std::endl; 
  std::cout << " 2 -- Sibson's C1 interpolant "<<std::endl;
  std::cout << " 3 -- Sibson's C1 interpolant without square-root"<<std::endl;
  std::cout << " 4 -- Farin's C1 interpolant" << std::endl;
  std::cin >> method;
  

  //INTERPOLATION:
  Coord_type value;
  std::pair<Coord_type, bool> res;
  std::vector< std::pair< Point_2, Coord_type > > coords;
  int n = points.size();
  ITraits traits;
  
  std::cout << "Interpolation at  "<<n  <<" grid points " << std::endl;
  for(int i=0;i<n;i++){
    Coord_type norm = 
      CGAL::natural_neighbor_coordinates_2(T, points[i],
					   std::back_inserter(coords)).second;
    assert(norm>0);  
    
    switch(method){
    case 0: 
      value = CGAL::linear_interpolation(coords.begin(),coords.end(),norm,
				       CGAL::Data_access<Point_value_map>(values)); 
      break;
    case 1: 
      res = CGAL::quadratic_interpolation(coords.begin(),coords.end(),
					  norm, points[i], 
					  CGAL::Data_access< Point_value_map>
					  (values),
					  CGAL::Data_access< Point_vector_map>
					  (gradients),traits); break;
    case 2:  
      res = CGAL::sibson_c1_interpolation(coords.begin(),coords.end(),
       					  norm, points[i], 
       					  CGAL::Data_access<Point_value_map>
       					  (values),
       					  CGAL::Data_access<Point_vector_map>
       					  (gradients), traits); break;
    case 3:  
      res = CGAL::sibson_c1_interpolation_square(coords.begin(),coords.end(),
						 norm, points[i], 
						 CGAL::Data_access<Point_value_map>
						 (values),
						 CGAL::Data_access<Point_vector_map>
						 (gradients), traits); break;
    case 4:  
      res = CGAL::farin_c1_interpolation(coords.begin(),coords.end(),
					 norm, points[i], 
					 CGAL::Data_access<Point_value_map>
					 (values),
					 CGAL::Data_access<Point_vector_map>
					 (gradients), traits); break;
      
    default: std::cout <<"No valid choice of interpolant." <<
		  std::endl; break;
    }
    if(method==0)
      points_3.push_back(Point_3(points[i].x(), points[i].y(),value));
    else if(res.second)
      points_3.push_back(Point_3(points[i].x(), points[i].y(),res.first));
    else std::cout <<"Interpolation failed"<<std::endl;
    coords.clear();
  }
  /************** end of Coordinate computation **************/
  //dump_off_file_quadrilateral(points_3, m);
  char ch;
  
  //viewer
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();
  gv.set_line_width(2);
  
  visu_graph(gv, points_3,m);

  std::cout << "The data points are displayed in blue in the geomview"
	    << " application." << std::endl;
  gv << CGAL::BLUE;
  visu_points(gv,sample_3);
  
  //show the gradients
  if(method>0){
    std::cout << "The function gradients are displayed by red lines "
	      <<" in the geomview application." << std::endl;
    gv <<CGAL::RED;
    gv << Segment_3(Point_3(h/3,h/3,w),Point_3(h/3,h/3,w)+ Vector_3(-g,-g,0));
    gv << Segment_3(Point_3(-h/3,h/3,w),Point_3(-h/3,h/3,w)+Vector_3(g,-g,0));
    gv << Segment_3(Point_3(-h/3,-h/3,w),Point_3(-h/3,-h/3,w)+Vector_3(g,g,0));
    gv << Segment_3(Point_3(h/3,-h/3,w),Point_3(h/3,-h/3,w)+Vector_3(-g,g,0));
  }
  switch(method){
  case 0: std::cout << " linear interpolation"; break;
  case 1: std::cout << " simple quadratic interpolant"; break;
  case 2: std::cout << " Sibson's C1 interpolant";break;
  case 3: std::cout << " Sibson's C1 interpolant without square-root"; break;
  case 4: std::cout << " Farin's C1 interpolant"; break;
  }
  std::cout << std::endl;
  std::cout << "Enter any character to quit." << std::endl;
  std::cin >> ch;
  
  return 0; 
}
#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
 
