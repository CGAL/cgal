
#define CGAL_FORCE_IFORMAT_DOUBLE

//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  Remeshing_triangulation tr;

  CGAL::Timer timer;
  timer.start();
  std::ifstream in("cheese.mesh");
  if (CGAL::SMDS_3::build_triangulation_from_file(in, tr))
    std::cout << "build triangulation ok" << std::endl;


  std::cout << "Read " << tr.number_of_vertices() << " vertices in " << timer.time() << " seconds.\n";

  return EXIT_SUCCESS;
}
