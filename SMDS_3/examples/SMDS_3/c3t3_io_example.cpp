#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/tags.h>

#include <CGAL/IO/C3t3_io.h>
#include <CGAL/IO/File_medit.h>

#include <fstream>


using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Subdomain_index = int;
using Surface_patch_index = unsigned char;
using Curve_index = char;
using Corner_index = short;

using Cb = CGAL::Simplicial_mesh_cell_base_3<K, Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                                  Curve_index, Corner_index>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Triangulation = CGAL::Triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;


int main(int argc, char* argv[])
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::string medit_elephant = (argc > 1)
                             ? std::string(argv[1])
                             : CGAL::data_file_path("meshes/elephant.mesh");

  Triangulation tr;
  std::ifstream is(medit_elephant, std::ios_base::in);
  if(!CGAL::IO::read_MEDIT(is, tr)) {
    std::cerr << "Failed to read " << medit_elephant << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = (argc > 1) ? std::string(argv[1])
                       : CGAL::data_file_path("meshes/elephant.c3t3.cgal");

  C3t3 c3t3;
  c3t3.set_triangulation(tr);

  const bool binary = false;
  std::ofstream os("after_load.c3t3.cgal");
  CGAL::IO::save_c3t3(os, c3t3, false);
  os.close();


  std::ifstream is2("after_load.c3t3.cgal");
  if(!CGAL::IO::load_c3t3(is2, c3t3, false))
  {
    std::cerr << "Failed to read" << std::endl;
    return EXIT_FAILURE;
  }

  // [do something]


  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
