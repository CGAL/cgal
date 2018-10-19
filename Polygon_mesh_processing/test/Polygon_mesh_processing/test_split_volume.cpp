#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <boost/core/ref.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if (argc==1)
  {
    std::cerr << "Nothing tested\n";
    return 1;
  }

  for(int i=1;i<argc; ++i)
  {
    std::cout << "Handling " << argv[i] << "\n";

    std::ifstream input(argv[i]);
    assert(input);
    Surface_mesh sm;
    input >> sm;

    std::vector<PMP::Volume_error_code> error_codes;
    std::vector< std::vector<std::size_t> > nested_cc_per_cc;

    Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> fccmap =
      sm.add_property_map<Surface_mesh::Face_index, std::size_t>("f:CC").first;

    Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> vol_id_map =
      sm.add_property_map<Surface_mesh::Face_index, std::size_t>().first;
    std::size_t nb_vol = PMP::volume_connected_components(
      sm, vol_id_map,
      CGAL::parameters::error_codes(boost::ref(error_codes)).
      face_connected_component_map(fccmap).
      volume_inclusions(boost::ref(nested_cc_per_cc)));

    std::cout << "  found " << nb_vol << " volumes\n";
    typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
    Filtered_graph vol_mesh(sm, 0, vol_id_map);
    for(std::size_t id = 0; id < nb_vol; ++id)
    {
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
}
