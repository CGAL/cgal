// to avoid ambiguities
#define LEDA_NO_MIN_MAX_TEMPL
// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/basic.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

typedef CGAL::leda_rat_kernel_traits           K;

typedef K::Point_3    Point;
typedef K::Segment_3  Segment;
typedef K::Ray_3      Ray;
typedef K::Line_3     Line;
typedef K::Triangle_3 Triangle;

typedef CGAL::leda_rat_kernel_traits            Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;

typedef CGAL::Triangulation_cell_base_3<Gt> Df;
typedef CGAL::Alpha_shape_cell_base_3<Gt, Df>  Fb;

typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Triangulation_3;

typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

typedef Alpha_shape_3::Cell  Cell;
typedef Alpha_shape_3::Vertex Vertex;
typedef Alpha_shape_3::Edge Edge;
typedef Alpha_shape_3::Cell_handle  Cell_handle;
typedef Alpha_shape_3::Vertex_handle Vertex_handle;

typedef Alpha_shape_3::Cell_circulator  Cell_circulator;

typedef Alpha_shape_3::Locate_type Locate_type;

typedef Alpha_shape_3::Cell_iterator  Cell_iterator;
typedef Alpha_shape_3::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_3::Edge_iterator  Edge_iterator;


typedef Alpha_shape_3::Coord_type Coord_type;
typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;

//---------------------------------------------------------------------

void
construct_alpha_shape(const std::list<Point> &V_p,
		      Alpha_shape_3::Mode mode,
		      Alpha_shape_3& A)
{ 
  std::vector<Gt::Segment_3> V_seg;
  
  int  n = A.make_alpha_shape(V_p.begin(), V_p.end());
  std::cout << "Inserted " << n  << " points" << std::endl;
  
  A.set_mode(mode);
}


void set_alpha(Alpha_shape_3& A, int alpha_index)
{
  // alpha corresponds to an index
  if (A.number_of_alphas() > 0)
  {
    if (alpha_index < 100)
      {
	int n = (alpha_index * A.number_of_alphas())/ 100;
	A.set_alpha(A.get_nth_alpha(n));
      }
    else
      {
	Alpha_shape_3::Alpha_iterator alpha_end_it = A.alpha_end();
	A.set_alpha((*(--alpha_end_it))+1);
      }
  }
  else
    A.set_alpha(0);
}
    
//------------------ main -------------------------------------------

int main()
{

  Alpha_shape_3 A;

  std::list<Point> L;
  L.push_back(Point(0,0,10,1)); L.push_back(Point(100,100,100,1));
  L.push_back(Point(50,0,0,1)); L.push_back(Point(100,200,70,1));
  L.push_back(Point(10,10,22,1)); L.push_back(Point(66,88,11,1));


  construct_alpha_shape(L,Alpha_shape_3::GENERAL,A);

  std::cout << "Alpha Shape computed" << std::endl;

  int n(50);


  if (n == 0)
    A.set_alpha(*A.find_optimal_alpha(2));
  else
    set_alpha(A,n);      

  return 0;
}
