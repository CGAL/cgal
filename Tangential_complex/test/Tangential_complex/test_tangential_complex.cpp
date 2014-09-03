#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/point_generators_3.h>

#include <fstream>

int main()
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<3> >  Kernel;
  typedef Kernel::Point_d                         Point;
  const int INTRINSIC_DIMENSION = 2;

  const int NUM_POINTS = 50;
  CGAL::Random_points_on_sphere_3<Point> generator(3.0);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
    points.push_back(*generator++);

  CGAL::Tangential_complex<Kernel, INTRINSIC_DIMENSION> tc(
                                                 points.begin(), points.end());
  tc.compute_tangential_complex();
  
  std::stringstream output_filename;
  output_filename << "data/test_tc_" << INTRINSIC_DIMENSION << ".off";
  std::ofstream off_stream(output_filename.str());
  tc.export_to_off(off_stream);

  return 0;
}