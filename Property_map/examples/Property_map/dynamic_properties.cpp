
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <string>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

int main()
{
  Mesh mesh;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0), mesh);

  typedef boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<std::string> >::type VertexNameMap;
  VertexNameMap vnm  = get(CGAL::dynamic_vertex_property_t<std::string>(), mesh);
  put(vnm, *(vertices(mesh).first), "Paris");

  std::cout << get(vnm, *(vertices(mesh).first)) << std::endl;

  typedef boost::property_map<Mesh, CGAL::dynamic_halfedge_property_t<double> >::type TrafficDensityMap;
  TrafficDensityMap tdm = get(CGAL::dynamic_halfedge_property_t<double>(), mesh);
  put(tdm, *(halfedges(mesh).first), 0.7);

  std::cout << get(tdm, *(halfedges(mesh).first)) << std::endl;

  return 0;
}

