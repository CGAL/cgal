/***********************************************************************

Takes a list of points and returns a list of segments corresponding to
the Alpha Shape.

************************************************************************/

#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

typedef double coord_type;

typedef CGAL::Simple_cartesian<coord_type>  K;

typedef K::Point_3    Point;
typedef K::Segment_3  Segment;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_3<K> Gt;
typedef Gt::Point Wpoint;
typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;

typedef CGAL::Triangulation_cell_base_3<Gt> Df;
typedef CGAL::Alpha_shape_cell_base_3<Gt, Df>  Fb;

typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds> Triangulation_3;

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
construct_alpha_shape(const std::list<Wpoint> &V_p,
		      Alpha_shape_3::Mode mode,
		      Alpha_shape_3& A)
  // Generate Alpha Shape
{ 
  std::vector<Segment> V_seg;
  
  int  n = A.make_alpha_shape(V_p.begin(), V_p.end());
  std::cout << "Inserted " << n  << " points" << std::endl;
  
  A.set_mode(mode);
}

bool
file_input(std::list<Wpoint>& L)
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
      L.push_back(Wpoint(p,5));
    }
  std::cout << "Points inserted" << std::endl;
  return true;
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


  std::list<Wpoint> L;

  file_input(L);  
  construct_alpha_shape(L,Alpha_shape_3::GENERAL,A);

  std::cout << "Alpha Shape computed" << std::endl;

  int n(50);


  if (n == 0)
    A.set_alpha(*A.find_optimal_alpha(2));
  else
    set_alpha(A,n);      
  //std::cout << A;


  return 0;
}
