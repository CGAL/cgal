#define CGAL_PROFILE

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Static_filters.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Simple_cartesian<CGAL::MP_Float>   K4;
typedef CGAL::Simple_cartesian<double>   K0;
typedef CGAL::Filtered_kernel<K0>        K1;
typedef K1         K2; // Static_filters is now included in Filtered_kernel.
//typedef CGAL::Static_filters<K1>         K2;
typedef CGAL::Kernel_checker<K2, K4,
              CGAL::Cartesian_converter<K2, K4> >     K3;

typedef K3::Point_2    Point_2;
typedef K3::Point_3    Point_3;

// Random double almost in [0;1].
double my_rand()
{
  // Ensure 53 random bits, not 48.
  return drand48() + drand48()/1024;
}

// Random point in unit square.
Point_2 my_rand_p2()
{
  return Point_2(my_rand(), my_rand());
}

// Random point on unit circle.
Point_2 circle_rand_p2()
{
  double dx = my_rand();
  double dy = my_rand();
  double n = std::sqrt(dx*dx+dy*dy);
  return Point_2(dx/n, dy/n);
}

// Random point in unit cube.
Point_3 my_rand_p3()
{
  return Point_3(my_rand(), my_rand(), my_rand());
}

// Random point on unit sphere.
Point_3 sphere_rand_p3()
{
  double dx = my_rand();
  double dy = my_rand();
  double dz = my_rand();
  double n = std::sqrt(dx*dx + dy*dy + dz*dz);
  return Point_3(dx/n, dy/n, dz/n);
}

// Perturbation with given maximum relative epsilon.
void perturb(Point_2 &p, double rel_eps)
{
  p = Point_2(p.x()*(1+drand48()*rel_eps),
              p.y()*(1+drand48()*rel_eps));
}

// Perturbation with given maximum relative epsilon.
void perturb(Point_3 &p, double rel_eps)
{
  p = Point_3(p.x()*(1+drand48()*rel_eps),
              p.y()*(1+drand48()*rel_eps),
              p.z()*(1+drand48()*rel_eps));
}

void
test_orientation_2()
{
  // First test with random points.
  Point_2 p = my_rand_p2();
  Point_2 q = my_rand_p2();
  Point_2 r = my_rand_p2();

  CGAL::CGALi::orientation(p, q, r, K3());
  CGAL::CGALi::orientation(q, r, p, K3());
  CGAL::CGALi::orientation(r, q, p, K3());

  // Then with collinear points (up to roundoff errors).
  r = p+(p-q)*my_rand();

  CGAL::CGALi::orientation(p, q, r, K3());
  CGAL::CGALi::orientation(q, r, p, K3());
  CGAL::CGALi::orientation(r, q, p, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::CGALi::orientation(p, q, r, K3());
  CGAL::CGALi::orientation(q, r, p, K3());
  CGAL::CGALi::orientation(r, q, p, K3());
}

void
test_orientation_3()
{
  // First test with random points.
  Point_3 p = my_rand_p3();
  Point_3 q = my_rand_p3();
  Point_3 r = my_rand_p3();
  Point_3 s = my_rand_p3();

  CGAL::CGALi::orientation(p, q, r, s, K3());
  CGAL::CGALi::orientation(q, r, s, p, K3());
  CGAL::CGALi::orientation(r, s, p, q, K3());
  CGAL::CGALi::orientation(s, p, q, r, K3());

  // Then with coplanar points (up to roundoff errors).
  s = p + (p-q)*my_rand() + (p-r)*my_rand();

  CGAL::CGALi::orientation(p, q, r, s, K3());
  CGAL::CGALi::orientation(q, r, s, p, K3());
  CGAL::CGALi::orientation(r, s, p, q, K3());
  CGAL::CGALi::orientation(s, p, q, r, K3());

  // Then with some perturbation.
  perturb(s, 1.0/(1<<20)/(1<<20)); // 2^-40

  CGAL::CGALi::orientation(p, q, r, s, K3());
  CGAL::CGALi::orientation(q, r, s, p, K3());
  CGAL::CGALi::orientation(r, s, p, q, K3());
  CGAL::CGALi::orientation(s, p, q, r, K3());
}

void
test_side_of_oriented_circle_2()
{
  // First test with random points.
  Point_2 p = my_rand_p2();
  Point_2 q = my_rand_p2();
  Point_2 r = my_rand_p2();
  Point_2 s = my_rand_p2();

  CGAL::CGALi::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::CGALi::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::CGALi::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::CGALi::side_of_oriented_circle(s, p, q, r, K3());

  // Then with cocircular points (up to roundoff errors).
  p = circle_rand_p2();
  q = circle_rand_p2();
  r = circle_rand_p2();
  s = circle_rand_p2();

  CGAL::CGALi::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::CGALi::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::CGALi::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::CGALi::side_of_oriented_circle(s, p, q, r, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::CGALi::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::CGALi::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::CGALi::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::CGALi::side_of_oriented_circle(s, p, q, r, K3());
}

void
test_side_of_oriented_sphere_3()
{
  // First test with random points.
  Point_3 p = my_rand_p3();
  Point_3 q = my_rand_p3();
  Point_3 r = my_rand_p3();
  Point_3 s = my_rand_p3();
  Point_3 t = my_rand_p3();

  CGAL::CGALi::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::CGALi::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::CGALi::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::CGALi::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::CGALi::side_of_oriented_sphere(t, s, p, q, r, K3());

  // Then with cospherical points (up to roundoff errors).
  p = sphere_rand_p3();
  q = sphere_rand_p3();
  r = sphere_rand_p3();
  s = sphere_rand_p3();
  t = sphere_rand_p3();

  CGAL::CGALi::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::CGALi::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::CGALi::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::CGALi::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::CGALi::side_of_oriented_sphere(t, s, p, q, r, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::CGALi::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::CGALi::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::CGALi::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::CGALi::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::CGALi::side_of_oriented_sphere(t, s, p, q, r, K3());
}

void compute_epsilons()
{
  std::cerr.precision(20);
  K2::Orientation_3::compute_epsilon();
  K2::Side_of_oriented_sphere_3::compute_epsilon();
}

int main(int argc, char **argv)
{
  compute_epsilons();

  int loops = (argc < 2) ? 1000 : atoi(argv[1]);

  std::cout << "#loops = " << loops << " (can be changed on the command line)" << std::endl;

  std::cout << "Testing Orientation_2" << std::endl;
  for(int i=0; i<loops; ++i)
    test_orientation_2();

  std::cout << "Testing Orientation_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_orientation_3();

  std::cout << "Testing Side_of_oriented_circle_2" << std::endl;
  for(int i=0; i<loops; ++i)
    test_side_of_oriented_circle_2();

  std::cout << "Testing Side_of_oriented_sphere_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_side_of_oriented_sphere_3();

  std::cerr << "Total number of IA failures = "
            << CGAL::Interval_nt<false>::number_of_failures << std::endl;

  return 0;
}
