// file: examples/Convex_hull_3/dynamic_hull_3_ex.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/copy_n.h>

#include <list>

typedef CGAL::Simple_cartesian<double> SK;
typedef CGAL::Filtered_kernel<SK> FK;
struct K : public FK {};

typedef K::Point_3                                  Point_3;
typedef CGAL::Delaunay_triangulation_3<K>           Delaunay;
typedef Delaunay::Vertex_handle                     Vertex_handle;

int main()
{
  CGAL::Random_points_in_sphere_3<Point_3> gen(100.0);
  std::list<Point_3>   points;

  // generate 250 points randomly on a sphere of radius 100.0 
  // and insert them into the triangulation
  CGAL::copy_n(gen, 250, std::back_inserter(points) );
  Delaunay T;
  T.insert(points.begin(), points.end());

  std::list<Vertex_handle>  vertices;
  T.incident_vertices(T.infinite_vertex(), std::back_inserter(vertices));
  std::cout << "This convex hull of the 250 points has " 
            << vertices.size() << " points on it." << std::endl;

  // remove 25 of the input points 
  std::list<Vertex_handle>::iterator v_set_it = vertices.begin();
  for (int i = 0; i < 25; i++)
  {
     T.remove(*v_set_it);
     v_set_it++;
  }

  vertices.clear();
  T.incident_vertices(T.infinite_vertex(), std::back_inserter(vertices));
  std::cout << "After removal of 25 points, there are "
            << vertices.size() << " points on the convex hull." << std::endl;
  return 0;
}
