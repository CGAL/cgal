/***********************************************************************

Prend une liste de points et renvoie une liste de segments
correspondant a l'Alpha Shape.

************************************************************************/



#include <CGAL/basic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream.h>

#include <list>

#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

#include <CGAL/Alpha_shape_euclidean_traits_2.h>

#include <CGAL/Alpha_shape_2.h>



typedef double coord_type;

typedef CGAL::Cartesian<coord_type>  Rep;

typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Ray_2<Rep>  Ray;
typedef CGAL::Line_2<Rep>  Line;
typedef CGAL::Triangle_2<Rep>  Triangle;

typedef CGAL::Alpha_shape_euclidean_traits_2<Rep> Gt;

typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;
typedef CGAL::Alpha_shape_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Dtriangulation_2;

typedef CGAL::Alpha_shape_2<Gt,Tds>  Alpha_shape_2;

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
Construst_Alpha_shape(const std::list<Point> &V_p,
		      const Coord_type &Alpha,
		      bool mode)
  // Generate Alpha Shape
{ 
  std::vector<Gt::Segment> V_seg;
  Alpha_shape_2 A;
  
  int  n = A.make_Alpha_shape(V_p.begin(), V_p.end());
  cerr << "Inserted " << n  << " points" << endl;
  
  if (mode) 
    { A.set_mode(Alpha_shape_2::GENERAL); } 
  else 
    { A.set_mode(Alpha_shape_2::REGULARIZED); };
  A.set_alpha(Alpha);

  V_seg << A;
  
  return V_seg;
}

bool
file_input(std::list<Point>& L)
{

  ifstream is("./data/fin", ios::in, filebuf::openprot);

  if(is.fail())
    {
      std::cerr << "unable to open file for input" << std::endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cerr << "Reading " << n << " points" << std::endl;
  Point p;
  for( ; n>0 ; n--)
    {
      is >> p;
      L.push_back(p);
    }
  std::cerr << "Points inserted" << std::endl;
  return true;
}
    
//------------------ main -------------------------------------------

int main(int argc,  char* argv[])
{
 std::list<Point> L;
 file_input(L);
 std::vector<Gt::Segment> V =
   Construst_Alpha_shape(L,10000.0,Alpha_shape_2::GENERAL);
 std::cerr << "Alpha Shape computed" << std::endl;
}
