#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <iostream>
#include <fstream>

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
  const char* filename = (argc > 1) ? argv[1] : "./data/elephant.off";
  std::ifstream in(filename);
  in >> sm;
   //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = sm.add_property_map<vertex_descriptor,double>("v:distance",0).first;

  //pass in the idt object and its vertex_distance_map
  Heat_method_idt hm_idt(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm_idt.add_source(source);
  hm_idt.estimate_geodesic_distances(vertex_distance);

  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " from " << source << std::endl;
  }

  return 0;
}
