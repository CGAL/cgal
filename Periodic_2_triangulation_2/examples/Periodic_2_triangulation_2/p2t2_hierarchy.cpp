#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_triangulation_hierarchy_2.h>
#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/point_generators_2.h>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K>         Gt;

typedef CGAL::Periodic_2_triangulation_vertex_base_2<Gt>            Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>            Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2<Gt>              Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt, Tds>          Dt;
typedef CGAL::Periodic_2_triangulation_hierarchy_2<Dt>              Triangulation;

typedef Triangulation::Point                                        Point;
typedef Triangulation::Iso_rectangle                                Iso_rectangle;

typedef CGAL::Creator_uniform_2<double, Point>                      Creator;

int main( )
{
  std::cout << "insertion of 1000 random points" << std::endl;
  Triangulation t(Iso_rectangle(-1,-1, 1,1));
  CGAL::Random_points_in_square_2<Point, Creator> g(1.);
  std::copy_n(g, 1000, std::back_inserter(t));

  //verbose mode of is_valid ; shows the number of vertices at each  level
  std::cout << "The number of vertices at successive levels" << std::endl;
  assert(t.is_valid(true));

  return 0;
}
