#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <fstream>
#include <list>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Alpha_shape_vertex_base_3<K>               Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>                 Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>      Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds,CGAL::Fast_location>  Delaunay;
typedef CGAL::Alpha_shape_3<Delaunay>                    Alpha_shape_3;

typedef K::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;
typedef Alpha_shape_3::NT                           NT;

int main()
{
  Delaunay dt;
  std::ifstream is("./data/bunny_1000");
  int n;
  is >> n;
  Point p;
  std::cout << n << " points read" << std::endl;
  for( ; n>0 ; n--) {
    is >> p;
    dt.insert(p);
  }
  std::cout << "Delaunay computed." << std::endl;

  // compute alpha shape
  Alpha_shape_3 as(dt);
  std::cout << "Alpha shape computed in REGULARIZED mode by defaut."
	    << std::endl;

   // find optimal alpha values
  Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Smallest alpha value to get a solid through data points is "
	    << alpha_solid << std::endl;
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);
  return 0;
}
