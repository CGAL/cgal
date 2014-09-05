#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/point_generators_3.h>

#include <fstream>

int main()
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<3> >  Kernel;
  typedef Kernel::Point_d                         Point;
  const int INTRINSIC_DIMENSION = 2;

#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 500;
#endif
  CGAL::default_random = CGAL::Random(0); // NO RANDOM
  CGAL::Random_points_on_sphere_3<Point> generator(3.0);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    points.push_back(*generator++);
    //points.push_back(Point((double)(rand()%10000)/5000, (double)(rand()%10000)/5000, 0.)); // CJTODO : plane
  }

  CGAL::Tangential_complex<Kernel, INTRINSIC_DIMENSION> tc(
                                                 points.begin(), points.end());
  tc.compute_tangential_complex();
  
  std::stringstream output_filename;
  output_filename << "data/test_tc_" << INTRINSIC_DIMENSION << ".off";
  std::ofstream off_stream(output_filename.str());
  tc.export_to_off(off_stream);

  return 0;
}