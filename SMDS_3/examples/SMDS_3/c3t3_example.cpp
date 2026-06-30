#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Simplicial_mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/IO/File_medit.h>
#include <fstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Triangulation = CGAL::Simplicial_mesh_triangulation_3<K>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;

int main(int argc, char* argv[])
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::string filename = (argc > 1) ? std::string(argv[1])
                       : CGAL::data_file_path("meshes/elephant.mesh");

  Triangulation tr;

  std::ifstream is(filename, std::ios_base::in);
  if(!CGAL::IO::read_MEDIT(is, tr))
  {
    std::cerr << "Failed to read" << std::endl;
    return EXIT_FAILURE;
  }

  C3t3 c3t3;
  c3t3.triangulation() = tr;

  std::ofstream os("out.mesh");
  CGAL::IO::write_MEDIT(os, tr, CGAL::parameters::all_vertices(true));
  os.close();

  return EXIT_SUCCESS;
}
