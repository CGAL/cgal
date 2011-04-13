/***********************************************************************

Takes a list of points and returns a list of segments corresponding to the
Alpha shape.

************************************************************************/


#include <CGAL/Cartesian.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream> // The standard says sstream
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_euclidean_traits_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

// Choose an exact number type.
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
typedef leda_integer coord_type;
#elif defined CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
typedef CGAL::Gmpz coord_type;
#else
#  include <CGAL/MP_Float.h>
typedef CGAL::MP_Float coord_type;
#endif

typedef CGAL::Cartesian<coord_type>  K;

typedef K::Point_2  Point;
typedef K::Segment_2  Segment;
typedef K::Ray_2  Ray;
typedef K::Line_2  Line;
typedef K::Triangle_2  Triangle;

typedef CGAL::Alpha_shape_euclidean_traits_2<K> Gt;

typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;

typedef CGAL::Triangulation_face_base_2<Gt> Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>  Fb;

typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;

typedef Alpha_shape_2::Face  Face;
typedef Alpha_shape_2::Vertex Vertex;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Face_handle  Face_handle;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef Alpha_shape_2::Face_circulator  Face_circulator;
typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape_2::Locate_type Locate_type;

typedef Alpha_shape_2::Face_iterator  Face_iterator;
typedef Alpha_shape_2::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator  Edge_iterator;
typedef Alpha_shape_2::Edge_circulator  Edge_circulator;

typedef Alpha_shape_2::Coord_type Coord_type;
typedef Alpha_shape_2::Alpha_iterator Alpha_iterator;

//---------------------------------------------------------------------

std::vector<Gt::Segment>
construct_alpha_shape(const std::list<Point> &V_p,
		      const Coord_type &Alpha,
		      bool mode)
  // Generate Alpha Shape
{ 
  std::vector<Gt::Segment> V_seg;
  Alpha_shape_2 A;
  
  int  n = A.make_alpha_shape(V_p.begin(), V_p.end());
  std::cout << "Inserted " << n  << " points" << std::endl;
  
  if (mode) 
    { A.set_mode(Alpha_shape_2::GENERAL); } 
  else 
    { A.set_mode(Alpha_shape_2::REGULARIZED); };
  A.set_alpha(Alpha);

  //V_seg << A;
  //return V_seg; 

  return A.op_vect_seg(V_seg);
}

bool
file_input(std::list<Point>& L)
{

  std::ifstream is("./data/fin", std::ios::in);

  if(is.fail())
    {
      std::cerr << "unable to open file for input" << std::endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  Point p;
  for( ; n>0 ; n--)
    {
      is >> p;
      L.push_back(p);
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}
    
//------------------ main -------------------------------------------

int main()
{
  std::list<Point> L;
  file_input(L);
  std::vector<Gt::Segment> V =
    construct_alpha_shape(L,Coord_type(10000),Alpha_shape_2::GENERAL);
  std::cout << "Alpha Shape computed" << std::endl;
  return 0;
}
