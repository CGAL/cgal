#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

#include <fstream>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef CGAL::Delaunay_triangulation<K> DT;

void test(int dim)
{
  std::stringstream input_filename;
  input_filename << "data/points_" << dim << ".cin";
  std::ifstream in(input_filename.str());

  DT::Point p;
  std::vector<DT::Point> points;

  int dim_from_file;
  in >> dim_from_file;
  while(in >> p)
    points.push_back(p);

  // Build the Regular Triangulation
  DT dt(dim_from_file);
  dt.insert(points.begin(), points.end());
  CGAL_assertion(dt.is_valid(true));

  // Export
  std::stringstream output_filename;
  output_filename << "data/dt_dim" << dim << ".off";
  std::ofstream off_stream(output_filename.str());
  CGAL::export_triangulation_to_off(off_stream, dt);
}

int main()
{
  //test(2);
  //test(3);
  test(10);
  return 0;
}
