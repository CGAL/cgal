#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K>             Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
typedef CGAL::Triangulation_face_base_2<K>               Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>      Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>            Dt;
typedef CGAL::Triangulation_hierarchy_2<Dt>              Triangulation;
typedef Triangulation::Point                             Point;
typedef CGAL::Creator_uniform_2<double, Point>            Creator;

int main( )
{
  std::cout << "insertion of 1000 random points" << std::endl;
  Triangulation t;
  CGAL::Random_points_in_square_2<Point, Creator> g(1.);
  CGAL::cpp11::copy_n( g, 1000, std::back_inserter(t));

  //verbose mode of is_valid ; shows the number of vertices at each  level
  std::cout << "The number of vertices at successive levels" << std::endl;
  assert(t.is_valid(true));

  return 0;
}
