
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>

#include <string>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

int main()
{
  Mesh mesh;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0), mesh);

  typedef boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<std::string> >::type VertexNameMap;
  VertexNameMap vnm  = get(CGAL::dynamic_vertex_property_t<std::string>(), mesh, std::string("default"));
  put(vnm, *(vertices(mesh).first), "Paris");

  assert(get(vnm, *(vertices(mesh).first))=="Paris");
  assert(get(vnm, *(std::next(vertices(mesh).first)))=="default");

  std::cout << get(vnm, *(vertices(mesh).first)) << std::endl;
  std::cout << get(vnm, *(std::next(vertices(mesh).first))) << std::endl;

  typedef boost::property_map<Mesh, CGAL::dynamic_halfedge_property_t<double> >::type TrafficDensityMap;
  TrafficDensityMap tdm = get(CGAL::dynamic_halfedge_property_t<double>(), mesh, -1.);
  put(tdm, *(halfedges(mesh).first), 0.7);

  assert(get(tdm, *(halfedges(mesh).first))==0.7);
  assert(get(tdm, *(std::next(halfedges(mesh).first)))==-1.);

  std::cout << get(tdm, *(halfedges(mesh).first)) << std::endl;
  std::cout << get(tdm, *(std::next(halfedges(mesh).first))) << std::endl;

  return 0;
}

