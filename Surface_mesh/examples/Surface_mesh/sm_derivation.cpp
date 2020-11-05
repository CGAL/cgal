#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

namespace My {

  struct Mesh: public CGAL::Surface_mesh<Point_3> {
    typedef CGAL::Surface_mesh<Point_3> Base;
    std::string name;
  };

} // namespace My


#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

int main()
{
  My::Mesh mesh;
  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,1), mesh);
  typedef boost::graph_traits<My::Mesh>::vertex_descriptor vertex_descriptor;

  typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;
  Point_property_map ppm = get(CGAL::vertex_point, mesh);

  for(vertex_descriptor vd : vertices(mesh)){
    if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
      std::cout << vd << " at " << get(ppm, vd) << std::endl;
    }
  }
  std::cout << CGAL::Polygon_mesh_processing::bbox(mesh) << std::endl;

  return 0;
}





