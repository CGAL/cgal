#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/tetrahedron_soup_to_lcc.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <array>

#include <CGAL/IO/MEDIT.h>

using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3>;
using Point_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.mesh");

  std::vector<Point_3> points;
  std::vector<std::array<std::size_t,4>> tetras;
  std::vector<int> subdomains;
  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Cannot open " << filename << std::endl;
    return 1;
  }

  CGAL::IO::read_MEDIT(in, points, tetras, subdomains);

  //filter out subdomain 0
  std::vector<std::array<std::size_t,4>> filtered_tetras;
  for (std::size_t i=0; i<tetras.size(); ++i)
    if (subdomains[i]!=0)
      filtered_tetras.push_back(tetras[i]);

  LCC lcc;
  CGAL::tetrahedron_soup_to_lcc(points, filtered_tetras, lcc);

  CGAL::draw(lcc);
}
