#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/generators.h>


template <class Traits>
struct Mesh: public CGAL::Polyhedron_3<Traits> {
  std::string name;
};

#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS typename Traits
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME Mesh<Traits>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Polyhedron_3<Traits>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

int main()
{
  Mesh<K> mesh;
  CGAL::make_triangle(Point_3(0,0,0), Point_3(1,0,0), Point_3(1,1,1), mesh);
  typedef boost::graph_traits<Mesh<K>>::vertex_descriptor vertex_descriptor;

  typedef boost::property_map<Mesh<K>,CGAL::vertex_point_t>::type Point_property_map;
  Point_property_map ppm = get(CGAL::vertex_point, mesh);

  for(vertex_descriptor vd : vertices(mesh)){
    if (vd != boost::graph_traits<Mesh<K>>::null_vertex()){
      std::cout << get(ppm, vd) << std::endl;
    }
  }
  std::cout << CGAL::Polygon_mesh_processing::bbox(mesh) << std::endl;

  return 0;
}





