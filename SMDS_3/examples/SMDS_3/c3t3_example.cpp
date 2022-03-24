#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing.h>

#include <CGAL/tags.h>

#include <CGAL/IO/File_medit.h>
#include <fstream>


using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Subdomain_index = int;
using Surface_patch_index = unsigned char;
using Curve_index = char;
using Corner_index = short;

using Cb = CGAL::Simplicial_mesh_cell_base_3<Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                                  Curve_index, Corner_index>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Triangulation = CGAL::Triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.mesh";

  Triangulation tr;

  std::ifstream is(filename, std::ios_base::in);
  CGAL::IO::read_MEDIT(is, tr);

  // [call a remeshing algorithm]

  std::ofstream os("after_remeshing.mesh");
  CGAL::IO::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}
