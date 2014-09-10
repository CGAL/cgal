// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/Kernel_traits.h>

#include <fstream>
#include <math.h>

#include <boost/random/random_number_generator.hpp>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 5000;
#endif

template <typename Point>
std::vector<Point> generate_points_on_sphere(double radius)
{
  CGAL::Random_points_on_sphere_3<Point> generator(radius);
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
    points.push_back(*generator++);
  return points;
}

// a = big radius, b = small radius
template <typename Point>
std::vector<Point> generate_points_on_klein_bottle(double a, double b)
{
  typedef CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT u = rng.get_double(0, 6.2832);
    FT v = rng.get_double(0, 6.2832);
    points.push_back(Kernel().construct_point_d_object()(
      (a + b*cos(v))*cos(u), 
      (a + b*cos(v))*sin(u),
      b*sin(v)*cos(u/2),
      b*sin(v)*sin(u/2)));
  }
  return points;
}

int main()
{
  const int INTRINSIC_DIMENSION = 2;
  const int AMBIENT_DIMENSION = 4;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> > Kernel;
  typedef Kernel::Point_d                                        Point;
 
  //CGAL::default_random = CGAL::Random(0); // NO RANDOM

#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  tbb::task_scheduler_init init(1);
# else
  tbb::task_scheduler_init init(10);
# endif
#endif

  //std::vector<Point> points = generate_points_on_sphere<Point>(3.0);
  std::vector<Point> points = generate_points_on_klein_bottle<Point>(4., 3.);

  CGAL::Tangential_complex<
    Kernel, 
    INTRINSIC_DIMENSION, 
    CGAL::Parallel_tag> tc(points.begin(), points.end());

  tc.compute_tangential_complex();
  
  std::stringstream output_filename;
  output_filename << "data/test_tc_" << INTRINSIC_DIMENSION << 
    << "_in_R" << AMBIENT_DIMENSION << ".off";
  std::ofstream off_stream(output_filename.str());
  tc.export_to_off(off_stream);

  return 0;
}