#define PROFILING_MODE

#include <CGAL/basic.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>
#include <CGAL/Timer.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>           Traits;
typedef Traits::FT                                                              NT;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef CGAL::Hyperbolic_octagon_translation_matrix<NT>                         Octagon_matrix;
typedef Triangulation::Point                                                    Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_original_octagon                                        Side_of_original_octagon;
typedef Triangulation::Face_iterator                                            Face_iterator;

typedef CGAL::Cartesian<double>::Point_2                                        Point_double;
typedef CGAL::Creator_uniform_2<double, Point_double >                          Creator;

long calls_apply_identity(0);
long calls_apply_non_identity(0);
long calls_append_identity(0);
long calls_append_non_identity(0);

long calls_predicate_identity(0);
long calls_predicate_non_identity(0);
double time_predicate_identity(0);
double time_predicate_non_identity(0);
double time_remove_dp(0);

int main(int argc, char** argv)
{
  int N, iters;
  if(argc < 2)
  {
    std::cout << "usage: " << argv[0] << " [number_of_points_to_insert] [optional: number_of_iterations]" << std::endl;
    std::cout << "Defaulting to 100k points, 10 iterations..." << std::endl;
    iters = 10;
    N = 100000;
  } else {
    N = atoi(argv[1]);
    if (argc < 3)
      iters = 10;
    else
      iters = atoi(argv[2]);
  }

  Side_of_original_octagon pred;

  std::cout << "---- for best results, make sure that you have compiled me in Release mode ----" << std::endl;

  double extime = 0.0;

  for(int exec = 1; exec <= iters; ++exec)
  {
    std::vector<Point> pts;
    CGAL::Random_points_in_disc_2<Point_double, Creator> g(0.85);

    int cnt = 0;
    do
    {
      Point_double pd = *(++g);
      if(pred(pd) != CGAL::ON_UNBOUNDED_SIDE)
      {
        Point pt = Point(pd.x(), pd.y());
        pts.push_back(pt);
        ++cnt;
      }
    }
    while(cnt < N);

    std::cout << "iteration " << exec << ": inserting into triangulation (rational dummy points)... "; std::cout.flush();
    Triangulation tr;

    CGAL::Timer tt;
    tt.start();
    tr.insert(pts.begin(), pts.end());
    tt.stop();
    std::cout << "DONE! (# of vertices = " << tr.number_of_vertices() << ", time = " << tt.time() << " secs)" << std::endl;
    extime += tt.time();

    int bfc = 0;
    for(Face_iterator fit = tr.faces_begin(); fit != tr.faces_end(); fit++)
    {
      if(!(fit->translation(0).is_identity() &&
           fit->translation(1).is_identity() &&
           fit->translation(2).is_identity()))
      {
        ++bfc;
      }
    }

    Triangulation::size_type Nf = tr.number_of_faces();
    double perc = double(bfc)/double(Nf) * 100.0;
    std::cout << "Total number of faces      : " << Nf << std::endl;
    std::cout << "Faces crossing the boundary: " << bfc << std::endl;
    std::cout << "Percentage                 : " << perc << std::endl;

    std::cout << "Triangulation is valid: " << (tr.is_valid() ? "YES" : "NO") << std::endl;
  }

  double avgtime = extime / double(iters);
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Average execution time over " << iters << " iterations: " << avgtime << " secs" << std::endl << std::endl;


  std::cout << "Calls to append resulting in     identity: " << calls_append_identity << std::endl;
  std::cout << "Calls to append resulting in non-identity: " << calls_append_non_identity << std::endl;
  std::cout << "Percentage                               : " << (double(calls_append_non_identity)/double(calls_append_non_identity+calls_append_identity)*100.0) << std::endl << std::endl;
  std::cout << "Calls to apply  with             identity: " << calls_apply_identity << std::endl;
  std::cout << "Calls to apply  with         non-identity: " << calls_apply_non_identity << std::endl;
  std::cout << "Percentage                               : " << double(calls_apply_non_identity)/(double(calls_apply_non_identity+calls_apply_identity)*100.0) << std::endl << std::endl;

  std::cout << "Predicate calls with     identity translations: " << calls_predicate_identity << std::endl;
  std::cout << "Predicate calls with non-identity translations: " << calls_predicate_non_identity << std::endl;
  std::cout << "Percentage                               : " << double(calls_predicate_non_identity)/double(calls_predicate_non_identity+calls_predicate_identity)*100.0 << std::endl << std::endl;

  std::cout << "Time in predicates with     identity translations: " << time_predicate_identity << std::endl;
  std::cout << "Time in predicates with non-identity translations: " << time_predicate_non_identity << std::endl;
  std::cout << "Percentage                                  : " << double(time_predicate_non_identity)/double(time_predicate_non_identity+time_predicate_identity)*100.0 << std::endl << std::endl;

  std::cout << "Time to remove dummy points                 : " << time_remove_dp << std::endl;

  return EXIT_SUCCESS;
}
