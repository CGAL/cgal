#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/IO/File_medit.h>

#include <fstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.mesh";

  Remeshing_triangulation tr;

  std::ifstream is(filename, std::ios_base::in);
  CGAL::IO::read_MEDIT(is, tr);

  // [call to a remeshing algorithm]

  std::ofstream os("after_remeshing.mesh");
  CGAL::IO::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}
