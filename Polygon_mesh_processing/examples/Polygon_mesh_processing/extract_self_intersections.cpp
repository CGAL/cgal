#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Real_timer.h>
#include <CGAL/tags.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");
  int k_ring = (argc > 2) ? std::stoi(argv[2]) : 3;
  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::cout << "Using parallel mode? " << std::is_same<CGAL::Parallel_if_available_tag, CGAL::Parallel_tag>::value << std::endl;

  CGAL::Real_timer timer;
  timer.start();


  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  PMP::self_intersections<CGAL::Parallel_if_available_tag>(faces(mesh), mesh, std::back_inserter(intersected_tris));

  if (intersected_tris.empty()) {
      std::cout << "No self intersections" << std::endl;
      return 0;
  }
  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;

  typedef Mesh::Property_map<face_descriptor,bool> Selection_pmap;

  Selection_pmap selection_pmap = mesh.add_property_map<face_descriptor,bool>("f:selection",false).first;


  std::set<face_descriptor> selection;
  for(auto p : intersected_tris){
    selection.insert(p.first);
    selection.insert(p.second);
    put(selection_pmap, p.first, true);
    put(selection_pmap, p.second, true);
  }


  CGAL::expand_face_selection(selection, mesh, k_ring, selection_pmap, std::inserter(selection, selection.end()));

  CGAL::Face_filtered_graph<Mesh> selection_mesh(mesh, selection);

  Mesh out;

  CGAL::copy_face_graph(selection_mesh, out);


  std::cout << timer.time() << " sec." << std::endl;

  std::ofstream os("out.off");
  os.precision(17);
  os << out << std::endl;

  return EXIT_SUCCESS;
}
