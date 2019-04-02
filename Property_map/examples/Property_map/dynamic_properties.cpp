
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <string>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 MyPoint_3;
//typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::Polyhedron_3<K> Mesh;

using namespace CGAL;

int main()
{
  Mesh mesh;
  CGAL::make_triangle(MyPoint_3(0,0,0),MyPoint_3(1,0,0),MyPoint_3(1,1,0), mesh);

  typedef boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<bool> >::type VertexNameMap;
  VertexNameMap vnm  = get(CGAL::dynamic_vertex_property_t<bool>(), mesh);
  put(vnm, *(vertices(mesh).first), true);
  
  std::cout << get(vnm, *(vertices(mesh).first)) << std::endl;
  
  typedef boost::property_map<Mesh, CGAL::dynamic_halfedge_property_t<bool> >::type TrafficDensityMap;
  TrafficDensityMap tdm = get(CGAL::dynamic_halfedge_property_t<bool>(), mesh);
  put(tdm, *(halfedges(mesh).first), false);
  
  std::cout << get(tdm, *(halfedges(mesh).first)) << std::endl;

  return 0;
}

