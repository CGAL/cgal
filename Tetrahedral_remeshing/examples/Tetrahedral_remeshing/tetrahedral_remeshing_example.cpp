#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> T3;

bool load_binary_triangulation(std::istream& is, T3& t3)
{
  std::string s;
  if (!(is >> s)) return false;
  bool binary = (s == "binary");
  if (binary) {
    if (!(is >> s)) return false;
  }
  if (s != "CGAL" || !(is >> s) || s != "c3t3")
    return false;

  std::getline(is, s);
  if (binary) CGAL::set_binary_mode(is);
  is >> t3;
  return bool(is);
}

bool save_binary_triangulation(std::ostream& os, const T3& t3)
{
  typedef T3::Geom_traits::FT FT;
  os << "binary CGAL c3t3\n";
  CGAL::set_binary_mode(os);
  return !!(os << t3);
}

int main(int argc, char* argv[])
{
  const char* filename     = (argc > 1) ? argv[1] : "data/triangulation_one_subdomain.binary.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in | std::ios::binary);

  T3 t3;
  if (!input)
    return false;

  if( !load_binary_triangulation(input, t3))
    return false;

  CGAL::tetrahedral_adaptive_remeshing(t3, target_edge_length);

  // save output
  const std::string file_in(filename);

  // binary
  std::string file_out = file_in.substr(0, file_in.find_first_of("."));
  file_out.append("_out.binary.cgal");
  std::ofstream out(file_out.c_str(), std::ios_base::out | std::ios_base::binary);
  save_binary_triangulation(out, t3);

  //// ascii
  //file_out = file_in.substr(0, file_in.find_first_of("."));
  //file_out.append("_out.mesh");
  //std::ofstream medit_out(file_out.c_str(), std::ios_base::out);
  //c3t3.output_to_medit(medit_out);

  return EXIT_SUCCESS;
}
