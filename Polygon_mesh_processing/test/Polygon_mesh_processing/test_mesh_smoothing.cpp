#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <boost/graph/graph_traits.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Kernel::Point_3 Point;



void test_angle_smoothing(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_angles(mesh);

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.7203429230262004);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 0);
      break;
    }
  }
}


void test_area_smoothing(const char* filename)
{
  std::cout.precision(17);
  std::ifstream input(filename);
  Mesh mesh;
  input >> mesh;
  input.close();

  boost::property_map<Mesh, CGAL::vertex_point_t>::type vpmap =
          get(CGAL::vertex_point, mesh);

  CGAL::Polygon_mesh_processing::smooth_areas(mesh);

  for (vertex_descriptor v : vertices(mesh))
  {
    if(!is_border(v, mesh))
    {
      Point p_c = get(vpmap, v);
      assert(p_c.x() == 0.6691415930575334);
      assert(p_c.y() == 0.5);
      assert(p_c.z() == 0);
      break;
    }
  }
}




int main(){

  const char* filename = "data/simple_polygon.off";
  test_angle_smoothing(filename);
  test_area_smoothing(filename);

  return 0;
}
