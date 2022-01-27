// Benchmark on the construction of Delaunay 2D with Epick (Kernel_23) vs
// something based on NewKernel_d that uses Eigen::Vector2d as Point_2.
#if  __cpp_aligned_new  >= 201606L
#if 1
#define CGAL_NEWKERNEL_D_USE_EIGEN_VECTOR 1
#include <CGAL/Epick_d.h>
#include <CGAL/NewKernel_d/Kernel_2_interface.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>

namespace CGAL {
struct Epick_2_help1
: Cartesian_filter_K<
    Cartesian_base_d<double, CGAL::Dimension_tag<2>>,
    Cartesian_base_d<Interval_nt_advanced, CGAL::Dimension_tag<2>>,
    Cartesian_base_d<internal::Exact_field_selector<double>::Type, CGAL::Dimension_tag<2>>
  >
{
};
struct Epick_2_help2
: Cartesian_filter_K<
    Epick_2_help1,
    Cartesian_base_d<Interval_nt_advanced, CGAL::Dimension_tag<2>>,
    Cartesian_base_d<internal::Exact_ring_selector<double>::Type, CGAL::Dimension_tag<2>>,
    Functors_without_division<CGAL::Dimension_tag<2>>::type
  >
{ };
struct Epick_2
: Kernel_2_interface<Cartesian_static_filters<CGAL::Dimension_tag<2>,Epick_2_help2>>
{ };
template <>
struct Triangulation_structural_filtering_traits<Epick_2> {
  typedef Tag_true Use_structural_filtering_tag;
};
}
typedef CGAL::Epick_2 K;
static_assert(std::is_same<K::Point_2, Eigen::Vector2d>::value, "");
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
#endif

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <random>
#include <stdlib.h>

typedef CGAL::Delaunay_triangulation_2<K>                    DT;
typedef DT::Point                                            Point_2;

int main(int argc,char** argv)
{
  int n = (argc > 1) ? atoi(argv[1]) : 100000;
  std::vector<Point_2> points;
  points.reserve( n );

  std::mt19937 gen(1234);
  std::uniform_real_distribution<> dis(1.0, 2.0);
  for(int i=0;i<n;++i)points.emplace_back(dis(gen),dis(gen));

  CGAL::Timer timer;
  timer.start();
  DT dt;
  dt.insert(points.begin(), points.end());
  timer.stop();
  std::size_t N = dt.number_of_vertices();

  std::cout << N << " points of type " << typeid(Point_2).name() << std::endl;
  std::cout << timer.time() << " seconds\n";

}
#else
#include <iostream>
int main(){
  std::cerr << "This program requires C++17 or later.\n";
}
#endif
