#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <boost/property_map/property_map.hpp>
#include <iostream>
#include <fstream>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef CGAL::Surface_mesh<Point>                           Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
typedef boost::graph_traits<Mesh>::faces_size_type          faces_size_type;

typedef Mesh::Property_map<face_descriptor, faces_size_type> FCCmap;
  typedef CGAL::Face_filtered_graph<Mesh> Filtered_graph;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  std::ifstream input((argc > 1) ? argv[1] : "data/blobby_3cc.off");

  Mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  FCCmap fccmap = mesh.add_property_map<face_descriptor, faces_size_type>("f:CC").first;

  faces_size_type num = PMP::connected_components(mesh,fccmap);

  std::cerr << "- The graph has " << num << " connected components (face connectivity)" << std::endl;

  Filtered_graph ffg(mesh, 0, fccmap);

  std::cout << "The faces in component 0 are:" << std::endl;
  for(boost::graph_traits<Filtered_graph>::face_descriptor f : faces(ffg)){
    std::cout  << f << std::endl;
  }

  if(num>1){
    std::vector<faces_size_type> components;
    components.push_back(0);
    components.push_back(1);

    ffg.set_selected_faces(components, fccmap);

    std::cout << "The faces in components 0 and 1 are:" << std::endl;
    for(Filtered_graph::face_descriptor f : faces(ffg)){
    std::cout  << f << std::endl;
    }
  }
  return 0;
}

