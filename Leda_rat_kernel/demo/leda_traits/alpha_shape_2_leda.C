// see example from CGAL Alpha shapes ...

#define LEDA_NO_MIN_MAX_TEMPL

#include <CGAL/basic.h>
#include <LEDA/rational.h>

// to remove ambiguity of min and max (std and leda)...
/*
namespace CGAL {

inline const leda_rational& min(const leda_rational& a, const leda_rational& b)
{ return (a<b) ? a : b; }

inline const leda_rational& max(const leda_rational& a, const leda_rational& b)
{ return (a>b) ? a : b; }

}
*/

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <iostream>
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_euclidean_traits_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

 
typedef CGAL::leda_rat_kernel_traits                K;
typedef K::Point_2                                  Point;
typedef K::Segment_2                                Segment;
typedef CGAL::leda_rat_kernel_traits                Gt;
typedef CGAL::Alpha_shape_vertex_base_2<Gt>         Vb;

typedef CGAL::Triangulation_face_base_2<Gt>         Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>       Fb;

typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>                 Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                   Alpha_shape_2;

typedef Alpha_shape_2::Vertex        Vertex;
typedef Alpha_shape_2::Edge          Edge;
typedef Alpha_shape_2::Face_handle   Face_handle;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;
typedef Alpha_shape_2::Locate_type        Locate_type;

typedef Alpha_shape_2::Vertex_iterator    Vertex_iterator;

typedef Alpha_shape_2::Coord_type         Coord_type;
typedef Alpha_shape_2::Alpha_iterator     Alpha_iterator;

//---------------------------------------------------------------------

std::vector<Gt::Segment>
construct_alpha_shape(const std::list<leda_rat_point>& V_p,
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
  return A.op_vect_seg(V_seg);
}


int main()
{
  std::list<leda_rat_point> L;

  //fill list with a few points ...
  L.push_back(Point(0,0,1)); L.push_back(Point(100,100,1));
  L.push_back(Point(50,0,1)); L.push_back(Point(100,200,1));
  L.push_back(Point(10,10,1)); L.push_back(Point(66,88,1));

  std::vector<Gt::Segment> V =
    construct_alpha_shape(L,Coord_type(10000),Alpha_shape_2::GENERAL);
  std::cout << "Alpha Shape computed" << std::endl;
  return 0;
}
