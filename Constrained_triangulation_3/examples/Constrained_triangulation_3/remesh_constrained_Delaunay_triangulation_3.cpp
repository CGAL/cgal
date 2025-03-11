#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE 1
#define CGAL_TETRAHEDRAL_REMESHING_DEBUG 1

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_triangulation_3.h>

#include <fstream>
#include <string>

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb  = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K>;
using Cb  = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Tr  = CGAL::Triangulation_3<K, Tds>;
using CCDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;

int main(int argc, char* argv[])
{
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mpi.off");
  std::ifstream in(filename);

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!in || !CGAL::IO::read_OFF(in, mesh)) {
    std::cerr << "Error: cannot read file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  CCDT ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CCDT>(mesh);
  //! [move ccdt to tr]
  Tr tr = std::move(ccdt).triangulation();
  //! [move ccdt to tr]
  std::cout << "Number of vertices in tr: "
            << tr.number_of_vertices() << std::endl;
  CGAL::tetrahedral_isotropic_remeshing(tr, 0.1,
    CGAL::parameters::number_of_iterations(3));

  std::cout << "Number of vertices in tr: "
            << tr.number_of_vertices() << std::endl;
  CGAL::draw(ccdt.triangulation());

  return EXIT_SUCCESS;
}
