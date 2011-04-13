//
// file: examples/Convex_hull_3/dynamic_hull_3_ex.C
//
#ifdef _MSC_VER
#define Cartesian Ca
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/copy_n.h>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
// NOTE: the choice of double here for a number type may cause problems
//       for degenerate point sets
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

#include <set>

typedef CGAL::Cartesian<RT>                         K;
typedef K::Point_3                                  Point_3;
typedef CGAL::Triangulation_cell_base_3<K>          Cb;
typedef CGAL::Triangulation_vertex_base_3<K>        Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>      Delaunay;
typedef Delaunay::Vertex_handle                     Vertex_handle;

int main()
{
  CGAL::Random_points_in_sphere_3<Point_3> gen(100.0);

  // generate 250 points randomly on a sphere of radius 100.0 
  // and insert them into the triangulation
  Delaunay T;
  CGAL::copy_n( gen, 250, std::back_inserter(T) );

  std::set<Vertex_handle>  vertices;
  T.incident_vertices(T.infinite_vertex(), vertices);
  std::cout << "This convex hull of the 250 points has " 
            << vertices.size() << " points on it." << std::endl;

  // remove 25 of the input points 
  std::set<Vertex_handle>::iterator v_set_it = vertices.begin();
  for (int i = 0; i < 25; i++)
  {
     T.remove(*v_set_it);
     v_set_it++;
  }

  vertices.clear();
  T.incident_vertices(T.infinite_vertex(), vertices);
  std::cout << "After removal of 25 points, there are "
            << vertices.size() << " points on the convex hull." << std::endl;
  return 0;
}
