// Checks that every internal variant of Aff_transformation_2/3 in the
// Cartesian kernels (identity, scaling, translation, rotation, reflection,
// general) behaves consistently, and that copies can be used concurrently.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

#include <cassert>
#include <thread>
#include <vector>
#include <iostream>

template <class K>
void test_2()
{
  typedef typename K::FT                   FT;
  typedef typename K::Point_2              Point_2;
  typedef typename K::Vector_2             Vector_2;
  typedef typename K::Direction_2          Direction_2;
  typedef typename K::Line_2               Line_2;
  typedef typename K::Circle_2             Circle_2;
  typedef typename K::Aff_transformation_2 T;

  std::vector<T> ts = {
    T(),
    T(CGAL::IDENTITY),
    T(CGAL::SCALING, FT(3)),
    T(CGAL::SCALING, FT(-2)),
    T(CGAL::TRANSLATION, Vector_2(1, 2)),
    T(CGAL::ROTATION, FT(3)/FT(5), FT(4)/FT(5)),
    T(CGAL::ROTATION, FT(3), FT(4), FT(5)),
    T(CGAL::REFLECTION, Line_2(Point_2(0, 1), Point_2(1, 1))),
    T(CGAL::REFLECTION, Line_2(Point_2(1, 0), Point_2(2, 3))),
    T(1, 2, 3, 4, 5, 6),
    T(1, 2, 3, 4, 5, 6, 7),
    T(0, -1, 1, 0),
    T(1, 1, 0, 0, 1, 0)                     // shear
  };

  const Point_2 p(1, -2);
  const Vector_2 v(2, 3);

  for(const T& a : ts)
  {
    // vectors are differences of points
    assert(a(p + v) - a(p) == a(v));
    // directions follow vectors
    assert(a(Direction_2(v)) == Direction_2(a(v)));
    // parity is the sign of the determinant of the linear part
    FT det = a.cartesian(0,0)*a.cartesian(1,1) - a.cartesian(0,1)*a.cartesian(1,0);
    assert(a.is_even() == CGAL::is_positive(det));
    assert(a.is_odd() == !a.is_even());
    // has_rotation must be true whenever the linear part is not diagonal
    if(a.cartesian(0,1) != 0 || a.cartesian(1,0) != 0)
      assert(a.has_rotation());
    // inverse
    assert(a.inverse()(a(p)) == p);
    // lines: the image of a line contains the images of its points
    Line_2 l(Point_2(0, 0), Point_2(1, 1));
    Line_2 la = a(l);
    assert(la.has_on(a(Point_2(0, 0))) && la.has_on(a(Point_2(1, 1))));
    // composition
    for(const T& b : ts)
      assert((a * b)(p) == a(b(p)));
  }

  // circle orientation: scaling by -1 is a rotation by 180 degrees in 2D
  Circle_2 c(Point_2(1, 0), FT(1), CGAL::COUNTERCLOCKWISE);
  assert(c.orthogonal_transform(T(CGAL::SCALING, FT(-1))).orientation() == CGAL::COUNTERCLOCKWISE);
  assert(c.orthogonal_transform(T(CGAL::REFLECTION, Line_2(Point_2(0,0), Point_2(1,0)))).orientation() == CGAL::CLOCKWISE);
}

template <class K>
void test_3()
{
  typedef typename K::FT                   FT;
  typedef typename K::Point_3              Point_3;
  typedef typename K::Vector_3             Vector_3;
  typedef typename K::Direction_3          Direction_3;
  typedef typename K::Plane_3              Plane_3;
  typedef typename K::Aff_transformation_3 T;

  std::vector<T> ts = {
    T(),
    T(CGAL::IDENTITY),
    T(CGAL::SCALING, FT(3)),
    T(CGAL::SCALING, FT(-2)),
    T(CGAL::SCALING, FT(6), FT(4)),
    T(CGAL::TRANSLATION, Vector_3(1, 2, 3)),
    T(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3),
    T(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 2),
    T(0, -1, 0, 1, 0, 0, 0, 0, 1),
    T(2, 0, 0, 0, 1, 0, 0, 0, -1, 3)
  };

  const Point_3 p(1, -2, 5);
  const Vector_3 v(2, 3, -1);

  for(const T& a : ts)
  {
    assert(a(p + v) - a(p) == a(v));
    assert(a(Direction_3(v)) == Direction_3(a(v)));
    FT det = CGAL::determinant(a.cartesian(0,0), a.cartesian(0,1), a.cartesian(0,2),
                               a.cartesian(1,0), a.cartesian(1,1), a.cartesian(1,2),
                               a.cartesian(2,0), a.cartesian(2,1), a.cartesian(2,2));
    assert(a.is_even() == CGAL::is_positive(det));
    if(a.cartesian(0,1) != 0 || a.cartesian(0,2) != 0 || a.cartesian(1,0) != 0 ||
       a.cartesian(1,2) != 0 || a.cartesian(2,0) != 0 || a.cartesian(2,1) != 0)
      assert(a.has_rotation());
    assert(a.inverse()(a(p)) == p);
    Plane_3 pl(1, 2, 3, 4);
    assert(a(pl).has_on(a(pl.point())));
    for(const T& b : ts)
      assert((a * b)(p) == a(b(p)));
  }
}

// Issue #9558: copies of a transformation used from several threads.
template <class K>
void test_threads()
{
  typedef typename K::Aff_transformation_3 T3;
  typedef typename K::Aff_transformation_2 T2;
  const T3 t3(CGAL::TRANSLATION, typename K::Vector_3(1, 2, 3));
  const T2 t2(CGAL::ROTATION, 3, 4, 5);
  std::vector<std::thread> threads;
  for(int i = 0; i < 8; ++i)
    threads.emplace_back([&]() {
      for(int j = 0; j < 10000; ++j) {
        T3 c3 = t3;   // copy
        T3 d3(c3);    // copy of copy, destroyed concurrently
        T2 c2 = t2;
        (void) d3(typename K::Point_3(j, 0, 0));
        (void) c2(typename K::Point_2(j, 0));
      }
    });
  for(std::thread& th : threads)
    th.join();
}

int main()
{
  typedef CGAL::Simple_cartesian<CGAL::Exact_rational> SC;
  typedef CGAL::Epeck                                   Epeck;

  test_2<SC>();
  test_3<SC>();
  test_2<Epeck>();
  test_3<Epeck>();
  test_threads<CGAL::Simple_cartesian<double> >();

  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}
