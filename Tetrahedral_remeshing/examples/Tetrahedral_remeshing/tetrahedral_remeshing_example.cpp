#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_binary_mesh_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> T3;
//todo : add specialization for Cell_base without info
// (does not compile with `void` instead of `int`)
typedef int Corner_index;
typedef int Curve_segment_index;

typedef CGAL::Mesh_complex_3_in_triangulation_3<T3, Corner_index, Curve_segment_index> C3t3;


int main(int argc, char* argv[])
{
  const char* filename     = (argc > 1) ? argv[1] : "data/triangulation_one_subdomain.binary.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in | std::ios::binary);

  C3t3 c3t3;
  if (!input)
    return false;

  if( !CGAL::Mesh_3::load_binary_file(input, c3t3))
    return false;

  CGAL::tetrahedral_adaptive_remeshing(c3t3.triangulation(), target_edge_length);

  // save output
  const std::string file_in(filename);

  // binary
  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append("_out.binary.cgal");
  std::ofstream out(file_out.c_str(), std::ios_base::out | std::ios_base::binary);
  CGAL::Mesh_3::save_binary_file(out, c3t3);

  // ascii
  file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append("_out.mesh");
  std::ofstream medit_out(file_out.c_str(), std::ios_base::out);
  c3t3.output_to_medit(medit_out);

  return EXIT_SUCCESS;
}
