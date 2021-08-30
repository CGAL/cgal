#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

#include <fstream>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef CGAL::Regular_triangulation<K> RT;

void test(int dim)
{
  std::stringstream input_filename;
  input_filename << "data/points_" << dim << ".cin";
  std::ifstream in(input_filename.str());

  RT::Weighted_point wp;
  std::vector<RT::Weighted_point> wpoints;

  int dim_from_file;
  in >> dim_from_file;
  while(in >> wp)
    wpoints.push_back(wp);

  // Build the Regular Triangulation
  RT rt(dim_from_file);
  rt.insert(wpoints.begin(), wpoints.end());
  CGAL_assertion(rt.is_valid(true));

  // Export
  std::stringstream output_filename;
  output_filename << "data/rt_dim" << dim << ".off";
  std::ofstream off_stream(output_filename.str());
  CGAL::export_triangulation_to_off(off_stream, rt);
}

int main()
{
  test(2);
  test(3);
  return 0;
}
