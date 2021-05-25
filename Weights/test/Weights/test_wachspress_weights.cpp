#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weights/wachspress_weights.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

int main() {

  // 2D configuration.
  const Point_2 t2 = Point_2(-1,  0);
  const Point_2 r2 = Point_2( 0, -1);
  const Point_2 p2 = Point_2( 1,  0);
  const Point_2 q2 = Point_2( 0,  0);

  // 3D configuration.
  const Point_3 t3 = Point_3(-1,  0, 1);
  const Point_3 r3 = Point_3( 0, -1, 1);
  const Point_3 p3 = Point_3( 1,  0, 1);
  const Point_3 q3 = Point_3( 0,  0, 1);

  // Compute weights.
  std::cout << "2D wachspress: " << CGAL::Weights::
    wachspress_weight(t2, r2, p2, q2) << std::endl;
  std::cout << "3D wachspress: " << CGAL::Weights::internal::
    wachspress_weight(t3, r3, p3, q3) << std::endl;

  // 2D configuration.
  const std::vector<Point_2> polygon2 = {t2, r2, p2, Point_2(0, 1)};

  std::vector<double> weights2;
  weights2.reserve(polygon2.size());
  CGAL::Weights::wachspress_weights_2(
    polygon2, q2, std::back_inserter(weights2));

  std::cout << "2D wachspress (polygon): ";
  for (const double weight2 : weights2)
    std::cout << weight2 << " ";
  std::cout << std::endl;

  // 3D configuration.
  CGAL::Triangulation_2_projection_traits_3<Kernel> ptraits(
    typename Kernel::Vector_3(0, 0, 1));

  const std::vector<Point_3> polygon3 = {t3, r3, p3, Point_3(0, 1, 1)};

  std::vector<double> weights3;
  weights3.reserve(polygon3.size());
  CGAL::Weights::wachspress_weights_2(
    polygon3, q3, std::back_inserter(weights3), ptraits);

  std::cout << "3D wachspress (polygon): ";
  for (const double weight3 : weights3)
    std::cout << weight3 << " ";
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

  // Merge it with the current test!

  // 2D configuration.
  // const Point_2 query2 = Point_2(+FT(0), FT(0));
  // const Point_2 vm2    = Point_2(+FT(1), FT(0));
  // const Point_2 vj2    = Point_2(+FT(0), FT(1));
  // const Point_2 vp2    = Point_2(-FT(1), FT(0));
  // Am = 0.5, Aj = 0.5, C = 1, B = 0.

  // 3D configuration 1. All points are coplanar on a horizontal plane.
  // const Point_3 query3 = Point_3(+FT(0), FT(0), FT(1));
  // const Point_3 vm3    = Point_3(+FT(1), FT(0), FT(1));
  // const Point_3 vj3    = Point_3(+FT(0), FT(1), FT(1));
  // const Point_3 vp3    = Point_3(-FT(1), FT(0), FT(1));
  // Am = 0.5, Aj = 0.5, C = 1, B = 0, 2D : 4.

  // 3D configuration 2. All points are coplanar on an arbitrary plane.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.812058795517);
  // const Point_3 vp3    = Point_3(+4.54574179746300, +7.9148892979530, +0.000000000000);
  // Am = 17.701, Aj = 26.574, C = 13.804, B = 30.471, 3D: 0.0293453.

  // 3D configuration 3. The first three points are coplanar whereas vp3 is slightly offsetted.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.0000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.0000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.8120587955170);
  // const Point_3 vp3    = Point_3(+4.22658807629900, +8.2522664687610, -0.2914612626125);

  // My area Am should coincide with the area below, but my areas Aj, C, and B should be
  // closer to the second configuration than areas below.
  // Am = 17.701, Aj = 26.641, C = 13.896, B = 30.524, 3D: 0.0293359.

  // 3D configuration 3. The first three points are coplanar whereas vp3 is even more slightly offsetted.
  // This configuration should converge to the second one comparing to the previous one.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.0000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.0000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.8120587955170);
  // const Point_3 vp3    = Point_3(+4.41756058270700, +8.0503895685670, -0.1170591355154);

  // My area Am should coincide with the area below, but my areas Aj, C, and B should be
  // closer to the second configuration than areas below.
  // Am = 17.701, Aj = 26.585, C = 13.819, B = 30.480, 3D: 0.0293438.

  // Compute weights.
  // WP wp;
  // std::cout << "2D wachspress: " << wp(query2, vm2, vj2, vp2) << std::endl;
  // std::cout << "3D wachspress: " << wp(query3, vm3, vj3, vp3) << std::endl;
