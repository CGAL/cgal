#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <fstream>
#include <iostream>

#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Polyhedron_3<Kernel>                           Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/elephant.off";
  std::ifstream in(filename);
  in >> sm;
  // map for the distance values to the source set
  boost::unordered_map<vertex_descriptor, double> vertex_distance;

  vertex_descriptor source = *(vertices(sm).first);
  
  CGAL::Heat_method_3::estimate_geodesic_distances(sm,
                                                   boost::make_assoc_property_map(vertex_distance),
                                                   source) ;

  std::cout << "Source vertex at: " << source->point() << std::endl;
  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd->point() << "  is at distance " << vertex_distance[vd] << std::endl;
  }

  return 0;
}
