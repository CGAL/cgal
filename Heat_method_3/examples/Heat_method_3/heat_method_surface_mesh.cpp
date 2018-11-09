#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <fstream>
#include <iostream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh> Heat_method;



int main(int argc, char* argv[])
{
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/sphere.off";
  std::ifstream in(filename);
  in >> sm;
  //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = sm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  Heat_method hm(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  hm.estimate_geodesic_distances(vertex_distance);
  
  Point_3 sp = sm.point(source);

  std::cout << "source: " << sp  << " " << source << std::endl;
  vertex_descriptor far;
  double sdistance = 0;
  
  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " from " << source << std::endl;
    if(squared_distance(sp,sm.point(vd)) > sdistance){
      far = vd;
      sdistance = squared_distance(sp,sm.point(vd));
    }
  }

  std::cout << "far: " << sm.point(far) << " " << far << std::endl;

  hm.add_source(far);
  hm.estimate_geodesic_distances(vertex_distance);

  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << "to the two sources" << std::endl;
  }

  std::cout << "done" << std::endl;
  return 0;
}
