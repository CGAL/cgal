#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/tags.h>

#include <iostream>
#include <fstream>

int main (int argc, char** argv)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;

  typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Open file
  std::ifstream in(argc > 1 ? argv[1] : "data/elephant.mesh",
                   std::ios_base::in);
  if(!in) {
    std::cerr << "Error! Cannot open file " << argv[1] << std::endl;
    return 1;
  }
  Tr tr;
  CGAL::IO::read_MEDIT(in, tr);

  std::ofstream os("elephant_out.mesh");
  CGAL::IO::write_MEDIT(os, tr);
  os.close();

  return EXIT_SUCCESS;
}
