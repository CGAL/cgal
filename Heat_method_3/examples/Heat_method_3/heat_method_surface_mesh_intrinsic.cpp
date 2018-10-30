#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>

#include <iostream>
#include <fstream>
#include <iostream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, CGAL::Heat_method_3::Intrinsic_Delaunay> Heat_method_idt;


int main(int argc, char* argv[])
{
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/bunny.off";
  std::ifstream in(filename);
  in >> sm;
  //the vertex distance map will hold the distance values from the source set
  Vertex_distance_map vdm_idt = sm.add_property_map<vertex_descriptor,double>("v:idt",0).first;

  //pass in the idt object and its vertex_distance_map
  Heat_method_idt hm_idt(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm_idt.add_source(source);
  hm_idt.fill_distance_map(vdm_idt);

  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(vdm_idt, vd) << " from " << source << std::endl;
  }

  return 0;
}
