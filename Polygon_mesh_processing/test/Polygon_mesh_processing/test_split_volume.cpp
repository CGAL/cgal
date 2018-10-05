#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

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

    Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> vol_id_map =
      sm.add_property_map<Surface_mesh::Face_index, std::size_t>().first;
    std::size_t nb_vol = 
      CGAL::Polygon_mesh_processing::volume_connected_components(sm, vol_id_map);
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
