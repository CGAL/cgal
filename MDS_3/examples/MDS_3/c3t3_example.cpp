#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/IO/File_medit.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.mesh";

  Remeshing_triangulation tr;

  std::ifstream is(filename, std::ios_base::in);
  CGAL::read_MEDIT(is, tr);

  // [call to a remeshing algorithm]

  std::ofstream os("after_remeshing.mesh");
  CGAL::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}
