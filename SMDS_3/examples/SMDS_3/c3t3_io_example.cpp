#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/tags.h>

#include <CGAL/IO/C3t3_io.h>
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

  std::string filename = (argc > 1) ? std::string(argv[1])
                       : CGAL::data_file_path("meshes/elephant.c3t3.cgal");

  C3t3 c3t3;

  std::ifstream is(filename, std::ios_base::binary);
  if(!CGAL::IO::load_c3t3(is, c3t3))
  {
    std::cerr << "Failed to read" << std::endl;
    return EXIT_FAILURE;
  }

  // [do something]

  std::ofstream os("after_load.c3t3.cgal", std::ios_base::binary);
  CGAL::IO::save_c3t3(os, c3t3);
  os.close();

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
