//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <iostream>
#include <fstream>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.mesh";
  const double target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1;

  Remeshing_triangulation tr;

  std::ifstream in(filename);
  if (CGAL::SMDS_3::build_triangulation_from_file(in, tr))
    std::cout << "build triangulation ok" << std::endl;
  else
    return EXIT_FAILURE;

  CGAL::Real_timer timer;
  timer.start();

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  timer.stop();
  std::cout << "Remeshing done in " << timer.time() << " seconds." << std::endl;

  // optional output, not written by the testsuite
  if(argc > 3)
  {
    std::ofstream out(argv[3]);
    CGAL::IO::write_MEDIT(out, tr);
  }

  return EXIT_SUCCESS;
}
