#define CGAL_LAZY_KERNEL_DEBUG
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/rational_rotation.h>
#include <cassert>
#include <thread>

int main()
{
  using K = CGAL::Exact_predicates_exact_constructions_kernel;
  using FT = K::FT;
  for(int sx : {-1, 1})
    for(int sy : {-1, 1})
    {
      FT sin, cos, den;
      const FT tolerance = FT(1) / FT(10000000000000000.0);
      CGAL::rational_rotation_approximation(FT(3 * sx), FT(4 * sy),
                                            sin, cos, den, tolerance, FT(1));
      assert(sin.depth() <= 1 && cos.depth() <= 1 && den.depth() <= 1);
      assert(CGAL::square(sin) + CGAL::square(cos) == CGAL::square(den));
      assert(CGAL::abs(sin / den - FT(4 * sy) / 5) < tolerance);
      assert(CGAL::abs(cos / den - FT(3 * sx) / 5) < tolerance);
    }

  std::thread worker([]()
  {
    const K::Point_2 a(-31420.8, -256500), b(-31419, -171900);
    const K::Vector_2 v(a, b);
    const K::Aff_transformation_2 rotation(CGAL::ROTATION,
                                         K::Direction_2(v.x(), -v.y()),
                                         FT(1), FT(1e16));
    const K::Vector_2 rotated = rotation.transform(v);
    assert(rotated.squared_length() == v.squared_length());
  });
  worker.join();
}
