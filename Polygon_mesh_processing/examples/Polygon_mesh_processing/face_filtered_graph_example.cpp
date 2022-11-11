#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <map>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;

typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor           face_descriptor;
typedef boost::graph_traits<Mesh>::faces_size_type           faces_size_type;

typedef Mesh::Property_map<face_descriptor, faces_size_type> FCCmap;
typedef CGAL::Face_filtered_graph<Mesh>                      Filtered_graph;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby_3cc.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  FCCmap fccmap = mesh.add_property_map<face_descriptor, faces_size_type>("f:CC").first;
  faces_size_type num = PMP::connected_components(mesh,fccmap);
  std::cerr << "- The graph has " << num << " connected components (face connectivity)" << std::endl;

  std::cout << "The faces in component 0 are:" << std::endl;
  Filtered_graph ffg(mesh, 0, fccmap);
  for(boost::graph_traits<Filtered_graph>::face_descriptor f : faces(ffg))
    std::cout << f << std::endl;

  if(num > 1)
  {
    std::vector<faces_size_type> components;
    components.push_back(0);
    components.push_back(1);

    std::cout << "The faces in components 0 and 1 are:" << std::endl;
    ffg.set_selected_faces(components, fccmap);
    for(Filtered_graph::face_descriptor f : faces(ffg))
      std::cout  << f << std::endl;
  }

  return 0;
}

