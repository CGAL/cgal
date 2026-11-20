#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/extrude_skeleton.h>

#include <type_traits>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using Float_kernel = CGAL::Simple_cartesian<float>;
using Double_kernel = CGAL::Simple_cartesian<double>;
using Long_double_kernel = CGAL::Simple_cartesian<long double>;

using EK = CGAL::Simple_cartesian<CGAL::Exact_rational>;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using EPECK_w_SQRT = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;

template <typename K>
using Kernel_selector = CGAL::Straight_skeleton_extrusion::internal::Extrusion_kernel_selector<K>;

static_assert(std::is_same_v<Kernel_selector<EPICK>::type, EPECK>);
static_assert(std::is_same_v<Kernel_selector<Float_kernel>::type, EPECK>);
static_assert(std::is_same_v<Kernel_selector<Double_kernel>::type, EPECK>);
static_assert(std::is_same_v<Kernel_selector<Long_double_kernel>::type, EPECK>);

static_assert(std::is_same_v<Kernel_selector<EK>::type, EK>);
static_assert(std::is_same_v<Kernel_selector<EPECK>::type, EPECK>);
static_assert(std::is_same_v<Kernel_selector<EPECK_w_SQRT>::type, EPECK_w_SQRT>);

int main(int, char**)
{
  std::cout << "OK" << std::endl;
  return EXIT_SUCCESS;
}
