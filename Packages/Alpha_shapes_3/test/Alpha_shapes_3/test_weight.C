/***********************************************************************

Takes a list of points and returns a list of segments corresponding to
the Alpha Shape.

************************************************************************/

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>


typedef double coord_type;
typedef CGAL::Simple_cartesian<coord_type>  SC;
typedef CGAL::Filtered_kernel<SC> K;

typedef K::Point_3    Point;
typedef K::Segment_3  Segment;
typedef K::Ray_3      Ray;
typedef K::Line_3     Line;
typedef K::Triangle_3 Triangle;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_3<K> Gt;
typedef Gt::Weighted_point   Weighted_point;

typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>   Fb;

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
typedef Alpha_shape_3::size_type   size_type;

typedef Alpha_shape_3::Finite_cells_iterator  Finite_cells_iterator;
typedef Alpha_shape_3::Finite_facets_iterator Finite_facets_iterator;
typedef Alpha_shape_3::Finite_edges_iterator  Finite_edges_iterator;
typedef Alpha_shape_3::Finite_vertices_iterator 
                                              Finite_vertices_iterator;
typedef Alpha_shape_3::Coord_type Coord_type;
typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;
typedef Alpha_shape_3::Alpha_shape_cells_iterator 
                                      Alpha_shape_cells_iterator;
typedef Alpha_shape_3::Alpha_shape_vertices_iterator 
                                      Alpha_shape_vertices_iterator;
typedef Alpha_shape_3::Alpha_shape_facets_iterator 
                                      Alpha_shape_facets_iterator;

//---------------------------------------------------------------------

#include <CGAL/count_alpha.C>

void
construct_alpha_shape(const std::list<Weighted_point> &V_p,
		      Alpha_shape_3::Mode mode,
		      Alpha_shape_3& A)
  // Generate Alpha Shape
{ 
  std::vector<Gt::Segment_3> V_seg;
  
  int  n = A.make_alpha_shape(V_p.begin(), V_p.end());
  std::cout << "Inserted " << n  << " points" << std::endl;
  
  A.set_mode(mode);
}

bool
file_input(std::list<Weighted_point>& L)
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
      L.push_back(Weighted_point(p,5*(n/10)));
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
  std::list<Weighted_point> L;
  bool verbose = false;
  
  // first a known small case
  // four groups of sphere :
  // - two groups of 4 intersecting spheres
  // - one group of three intersecting sphere
  // - one group of two intersecting sphere
  // the four groups are disjoint
  // Check specially the  $0$-shape which is the nerve of the union
  L.push_back(Weighted_point(Point(0.,0.,0.), 4));
  L.push_back(Weighted_point(Point(2.,2.,0.), 4));
  L.push_back(Weighted_point(Point(2.,0.,2.), 4));
  L.push_back(Weighted_point(Point(0.,2.,2.), 4));

  L.push_back(Weighted_point(Point(10.,0.,0.), 4));
  L.push_back(Weighted_point(Point(12.,2.,0.), 4));
  L.push_back(Weighted_point(Point(12.,0.,2.), 4));
  L.push_back(Weighted_point(Point(10.,2.,2.), 4));

  L.push_back(Weighted_point(Point(20.,0.,0.), 3));
  L.push_back(Weighted_point(Point(22.,0.,0.), 3));
  L.push_back(Weighted_point(Point(20.,2.,0.), 3));

  L.push_back(Weighted_point(Point(30.,0.,0.), 3));
  L.push_back(Weighted_point(Point(32.,0.,0.), 3));
  


  Alpha_shape_3 A( L.begin(), L.end(), 0, Alpha_shape_3::REGULARIZED);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);

  A.set_mode(Alpha_shape_3::GENERAL);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);

  assert(A.number_of_solid_components(0.) == 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(2))) <= 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);


  // test a bigger Alpha_shapes
  A.clear();
  L.clear();
  file_input(L);  
  construct_alpha_shape(L,Alpha_shape_3::GENERAL,A);
  std::cout << "Alpha Shape computed" << std::endl;
  std::cout << " test number_of_components - find_optimal_alpha "<<
    std::endl;
  A.set_alpha(*A.find_optimal_alpha(2));
  //show_alpha_values(A);
  // std::cerr << "optimal alpha " << *A.find_optimal_alpha(2) << std::endl;
  assert( A.number_of_solid_components() <= 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);
  return 0;
}
