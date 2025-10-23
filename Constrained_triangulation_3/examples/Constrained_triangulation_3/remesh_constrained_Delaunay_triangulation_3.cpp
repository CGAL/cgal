#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/property_map.h>

#include <unordered_set>
#include <fstream>
#include <string>

using K    = CGAL::Exact_predicates_inexact_constructions_kernel;

using Vbb  = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Vb   = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vbb>;

using Cbb  = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using Cb   = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cbb>;

using Tds  = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Tr   = CGAL::Triangulation_3<K, Tds>;
using CCDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;

// Triangulation for Remeshing
using CCDT_Tr = CCDT::Triangulation;
using Triangulation_3 = CGAL::Triangulation_3<K, CCDT_Tr::Triangulation_data_structure>;

using Vertex_handle = Triangulation_3::Vertex_handle;
using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;


int main(int argc, char* argv[])
{
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mpi.off");
  double target_edge_length = (argc > 2) ? std::stod(argv[2]) : 1.0;
  unsigned int iterations = (argc > 3) ? std::stoi(argv[3]) : 3;

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  CCDT ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CCDT>(mesh);

  std::ofstream out("ccdt.mesh");
  CGAL::IO::write_MEDIT(out, ccdt);
  out.close();

  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  namespace np = CGAL::parameters;
  namespace Tet_remesh = CGAL::Tetrahedral_remeshing;
  Tr tr = Tet_remesh::get_remeshing_triangulation(std::move(ccdt),
                                                  np::edge_is_constrained_map(constraints_pmap));
  std::cout << "There are " << tr.number_of_vertices() << " vertices in the constrained triangulation" << std::endl;

  CGAL::tetrahedral_isotropic_remeshing(tr,
                                        target_edge_length,
                                        np::number_of_iterations(iterations)
                                           .edge_is_constrained_map(constraints_pmap));

  std::cout << "There are " << tr.number_of_vertices() << " vertices after remeshing" << std::endl;

  std::ofstream ofs("remeshed.mesh");
  CGAL::IO::write_MEDIT(ofs, tr);

  return EXIT_SUCCESS;
}
