//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;


int main(int argc, char* argv[])
{
  const char* filename     = (argc > 1) ? argv[1] : "data/triangulation_one_subdomain.binary.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in | std::ios::binary);

  Remeshing_triangulation t3;
  if (!input)
    return EXIT_FAILURE;

  if( !CGAL::load_triangulation(input, t3))
    return EXIT_FAILURE;

  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length);

  // save output
  const std::string file_in(filename);
  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append("_out.binary.cgal");
  std::ofstream out(file_out.c_str(), std::ios_base::out | std::ios_base::binary);
  CGAL::save_binary_triangulation(out, t3);

  return EXIT_SUCCESS;
}
