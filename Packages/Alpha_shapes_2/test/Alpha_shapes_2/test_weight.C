/***********************************************************************

Prend une liste de points et renvoie une liste de segments
correspondant a l'Alpha Shape ponderee.

************************************************************************/



#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream>
// en fait le std dit sstream

#include <vector>
#include <list>

#include <CGAL/Weighted_point.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_2.h>

#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

//Choose the better number type as possible
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer coord_type;
#else//CGAL_USE_LEDA
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz coord_type;
#else//CGAL_USE_GMP
#include <CGAL/double.h>
typedef double coord_type;
#endif//CGAL_USE_GMP
#endif//CGAL_USE_LEDA

typedef CGAL::Cartesian<coord_type>  CRep;
typedef CGAL::Point_2<CRep> Point_base;
typedef CGAL::Weighted_point<Point_base,coord_type>  Point;
typedef CGAL::Segment_2<CRep>  Segment;
typedef CGAL::Ray_2<CRep>  Ray;
typedef CGAL::Line_2<CRep>  Line;
typedef CGAL::Triangle_2<CRep>  Triangle;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_2<CRep> Gt;

typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;

typedef CGAL::Regular_triangulation_face_base_2<Gt> Rf;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf>  Fb;

typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_2<Gt,Tds> Triangulation_2;

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
Construct_Alpha_shape(const std::list<Point> &V_p,
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

  //  V_seg << A;
  
  return V_seg;
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
      L.push_back(Point (p.point(),Coord_type(10)));
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}
    
//------------------ main -------------------------------------------

int main(int argc,  char* argv[])
{
  std::list<Point> L;
  file_input(L);
  std::vector<Gt::Segment> V =
    Construct_Alpha_shape(L,Coord_type(10000),Alpha_shape_2::GENERAL);
  std::cout << "Weighted Alpha Shape computed" << std::endl;
  return 0;
}
