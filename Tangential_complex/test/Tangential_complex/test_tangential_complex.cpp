//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

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
#include <CGAL/Mesh_3/Profiling_tools.h>

#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

#ifdef _DEBUG
  const int NUM_POINTS = 6;
#else
  const int NUM_POINTS = 50000;
#endif

template <typename Point>
std::vector<Point> generate_points_on_plane()
{
  typedef typename CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT x = rng.get_double(0, 5);
    FT y = rng.get_double(0, 5);
    points.push_back(Kernel().construct_point_d_object()(x, y, 0));
  }
  return points;
}

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
std::vector<Point> generate_points_on_klein_bottle_3D(
  double a, double b, bool uniform = false)
{
  typedef typename CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  // if uniform
  int num_lines = (int)sqrt(NUM_POINTS);
  int num_cols = NUM_POINTS/num_lines + 1;

  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      int k1 = i / num_lines;
      int k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    }
    else
    { 
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    double tmp = cos(u/2)*sin(v) - sin(u/2)*sin(2.*v);
    points.push_back(Kernel().construct_point_d_object()(
      (a + b*tmp)*cos(u), 
      (a + b*tmp)*sin(u),
      b*(sin(u/2)*sin(v) + cos(u/2)*sin(2.*v))));
  }
  return points;
}

// a = big radius, b = small radius
template <typename Point>
std::vector<Point> generate_points_on_klein_bottle_4D(
  double a, double b, bool uniform = false)
{
  typedef typename CGAL::Kernel_traits<Point>::type Kernel;
  typedef typename Kernel::FT FT;
  CGAL::Random rng;

  // if uniform
  int num_lines = (int)sqrt(NUM_POINTS);
  int num_cols = NUM_POINTS/num_lines + 1;

  std::vector<Point> points;
  points.reserve(NUM_POINTS);
  for (int i = 0 ; i != NUM_POINTS ; ++i)
  {
    FT u, v;
    if (uniform)
    {
      int k1 = i / num_lines;
      int k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    }
    else
    { 
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
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

  int i = 0;
  bool stop = false;
  //for ( ; !stop ; ++i)
  {
    Wall_clock_timer t;
    CGAL::default_random = CGAL::Random(i);
    std::cerr << "Random seed = " << i << std::endl;
  
    //std::vector<Point> points = generate_points_on_plane<Point>();
    //std::vector<Point> points = generate_points_on_sphere<Point>(3.0);
    //std::vector<Point> points = generate_points_on_klein_bottle_3D<Point>(4., 3.);
    std::vector<Point> points = generate_points_on_klein_bottle_4D<Point>(4., 3.);

    CGAL::Tangential_complex<
      Kernel, 
      INTRINSIC_DIMENSION, 
      CGAL::Parallel_tag> tc(points.begin(), points.end());
    double init_time = t.elapsed(); t.reset();

    tc.compute_tangential_complex();
    double computation_time = t.elapsed(); t.reset();

    std::set<std::set<std::size_t> > incorrect_simplices;
    //stop = !tc.check_if_all_simplices_are_in_the_ambient_delaunay(&incorrect_simplices);

    std::stringstream output_filename;
    output_filename << "data/test_tc_" << INTRINSIC_DIMENSION
      << "_in_R" << AMBIENT_DIMENSION << ".off";
    std::ofstream off_stream(output_filename.str());
    tc.export_to_off(off_stream, true, &incorrect_simplices, true);
    double export_time = t.elapsed(); t.reset();

    std::cerr << std::endl
      << "================================================" << std::endl
      << "Computation times (seconds): " << std::endl
      << "  * Tangential complex: " << init_time + computation_time
      << std::endl
      << "    - Init + kd-tree = " << init_time << std::endl
      << "    - TC computation = " << computation_time << std::endl
      << "  * Export to OFF: " << export_time << std::endl
      << "================================================" << std::endl
      << std::endl;
  }

  return 0;
}