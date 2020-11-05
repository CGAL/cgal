#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/blobby.off";
  std::ifstream input(filename);
  assert(input);
  Surface_mesh sm;
  input >> sm;

  // property map to assign a volume id for each face
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> vol_id_map =
    sm.add_property_map<Surface_mesh::Face_index, std::size_t>().first;

  std::vector<PMP::Volume_error_code> err_codes;
  // fill the volume id map
  std::size_t nb_vol =
    PMP::volume_connected_components(sm,
                                     vol_id_map,
                                     params::error_codes(std::ref(err_codes)));

  std::cout << "Found " << nb_vol << " volumes\n";

  // write each volume in an OFF file
  typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
  Filtered_graph vol_mesh(sm, 0, vol_id_map);
  for(std::size_t id = 0; id < nb_vol; ++id)
  {
    if (err_codes[id] != PMP::VALID_VOLUME)
      std::cerr << "There is an issue with volume #" << id << "\n";
    if(id > 0)
      vol_mesh.set_selected_faces(id, vol_id_map);
    Surface_mesh out;
    CGAL::copy_face_graph(vol_mesh, out);
    std::ostringstream oss;
    oss << "vol_" << id <<".off";
    std::ofstream os(oss.str().data());
    os << out;
  }
}
