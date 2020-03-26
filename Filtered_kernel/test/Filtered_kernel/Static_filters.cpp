#define CGAL_PROFILE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Random.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double>   K0;
typedef CGAL::Filtered_kernel<K0, true>  K2;
typedef CGAL::Filtered_kernel<K0, false> K4;
typedef CGAL::Kernel_checker<K2, K4>     K3;

typedef K3::Point_2    Point_2;
typedef K3::Point_3    Point_3;

CGAL::Random *r;

double rand_base()
{
  return r->get_double(0, 1);
}

// Random double almost in [0;1].
double my_rand()
{
  // Ensure 53 random bits, not 48.
  return rand_base() + rand_base()/1024;
}

// Random point in unit square.
Point_2 my_rand_p2()
{
  double x = my_rand(), y = my_rand();
  return Point_2(K2::Point_2(x, y), K4::Point_2(x, y));
}

// Random point on unit circle.
Point_2 circle_rand_p2()
{
  double dx = my_rand();
  double dy = my_rand();
  double n = std::sqrt(dx*dx+dy*dy);
  double x = dx/n, y = dy/n;
  return Point_2(K2::Point_2(x, y), K4::Point_2(x, y));
}

// Random point in unit cube.
Point_3 my_rand_p3()
{
  double x = my_rand(), y = my_rand(), z = my_rand();
  return Point_3(K2::Point_3(x, y, z), K4::Point_3(x, y, z));
}

// Random point on unit sphere.
Point_3 sphere_rand_p3()
{
  double dx = my_rand();
  double dy = my_rand();
  double dz = my_rand();
  double n = std::sqrt(dx*dx + dy*dy + dz*dz);
  double x = dx/n, y = dy/n, z = dz/n;
  return Point_3(K2::Point_3(x, y, z), K4::Point_3(x, y, z));
}

// Perturbation with given maximum relative epsilon.
void perturb(Point_2 &p, double rel_eps)
{
  double x = p.first.x()*(1+rand_base()*rel_eps);
  double y = p.first.y()*(1+rand_base()*rel_eps);
  p = Point_2(K2::Point_2(x, y), K4::Point_2(x, y));
}

// Perturbation with given maximum relative epsilon.
void perturb(Point_3 &p, double rel_eps)
{
  double x = p.first.x()*(1+rand_base()*rel_eps);
  double y = p.first.y()*(1+rand_base()*rel_eps);
  double z = p.first.z()*(1+rand_base()*rel_eps);
  p = Point_3(K2::Point_3(x, y, z), K4::Point_3(x, y, z));
}

// Pick a random point on the segment [p, q].
Point_2 pick_collinear(const Point_2 &p, const Point_2 &q)
{
  // Kernel_checker currently prevents from using the following :
  // r = p+(p-q)*my_rand();
  double r1 = my_rand();
  double x = p.first.x() + (q.first.x() - p.first.x())*r1;
  double y = p.first.y() + (q.first.y() - p.first.y())*r1;
  return Point_2(K2::Point_2(x, y), K4::Point_2(x, y));
}

// Pick a random point on the triangle [p, q, r].
Point_3 pick_coplanar(const Point_3 &p, const Point_3 &q, const Point_3 &r)
{
  // Kernel_checker currently prevents from using the following :
  // s = p + (p-q)*my_rand() + (p-r)*my_rand();  (almost)
  double r1 = my_rand(), r2 = my_rand();
  double x = p.first.x() + (q.first.x() - p.first.x())*r1
                         + (r.first.x() - p.first.x())*r2;
  double y = p.first.y() + (q.first.y() - p.first.y())*r1
                         + (r.first.y() - p.first.y())*r2;
  double z = p.first.z() + (q.first.z() - p.first.z())*r1
                         + (r.first.z() - p.first.z())*r2;
  return Point_3(K2::Point_3(x, y, z), K4::Point_3(x, y, z));
}

void
test_orientation_2()
{
  // First test with random points.
  Point_2 p = my_rand_p2();
  Point_2 q = my_rand_p2();
  Point_2 r = my_rand_p2();

  CGAL::internal::orientation(p, q, r, K3());
  CGAL::internal::orientation(q, r, p, K3());
  CGAL::internal::orientation(r, q, p, K3());

  // Then with collinear points (up to roundoff errors).
  // r = p+(p-q)*my_rand();
  r = pick_collinear(p, q);

  CGAL::internal::orientation(p, q, r, K3());
  CGAL::internal::orientation(q, r, p, K3());
  CGAL::internal::orientation(r, q, p, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::internal::orientation(p, q, r, K3());
  CGAL::internal::orientation(q, r, p, K3());
  CGAL::internal::orientation(r, q, p, K3());
}

void
test_orientation_3()
{
  // First test with random points.
  Point_3 p = my_rand_p3();
  Point_3 q = my_rand_p3();
  Point_3 r = my_rand_p3();
  Point_3 s = my_rand_p3();

  CGAL::internal::orientation(p, q, r, s, K3());
  CGAL::internal::orientation(q, r, s, p, K3());
  CGAL::internal::orientation(r, s, p, q, K3());
  CGAL::internal::orientation(s, p, q, r, K3());

  // Then with coplanar points (up to roundoff errors).
  //s = p + (p-q)*my_rand() + (p-r)*my_rand();
  s = pick_coplanar(p, q, r);

  CGAL::internal::orientation(p, q, r, s, K3());
  CGAL::internal::orientation(q, r, s, p, K3());
  CGAL::internal::orientation(r, s, p, q, K3());
  CGAL::internal::orientation(s, p, q, r, K3());

  // Then with some perturbation.
  perturb(s, 1.0/(1<<20)/(1<<20)); // 2^-40

  CGAL::internal::orientation(p, q, r, s, K3());
  CGAL::internal::orientation(q, r, s, p, K3());
  CGAL::internal::orientation(r, s, p, q, K3());
  CGAL::internal::orientation(s, p, q, r, K3());
}

void
test_side_of_oriented_circle_2()
{
  // First test with random points.
  Point_2 p = my_rand_p2();
  Point_2 q = my_rand_p2();
  Point_2 r = my_rand_p2();
  Point_2 s = my_rand_p2();

  CGAL::internal::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::internal::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::internal::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::internal::side_of_oriented_circle(s, p, q, r, K3());

  // Then with cocircular points (up to roundoff errors).
  p = circle_rand_p2();
  q = circle_rand_p2();
  r = circle_rand_p2();
  s = circle_rand_p2();

  CGAL::internal::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::internal::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::internal::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::internal::side_of_oriented_circle(s, p, q, r, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::internal::side_of_oriented_circle(p, q, r, s, K3());
  CGAL::internal::side_of_oriented_circle(q, r, s, p, K3());
  CGAL::internal::side_of_oriented_circle(r, s, p, q, K3());
  CGAL::internal::side_of_oriented_circle(s, p, q, r, K3());
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

  CGAL::internal::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::internal::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::internal::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::internal::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::internal::side_of_oriented_sphere(t, s, p, q, r, K3());

  // Then with cospherical points (up to roundoff errors).
  p = sphere_rand_p3();
  q = sphere_rand_p3();
  r = sphere_rand_p3();
  s = sphere_rand_p3();
  t = sphere_rand_p3();

  CGAL::internal::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::internal::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::internal::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::internal::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::internal::side_of_oriented_sphere(t, s, p, q, r, K3());

  // Then with some perturbation.
  perturb(r, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::internal::side_of_oriented_sphere(p, q, r, s, t, K3());
  CGAL::internal::side_of_oriented_sphere(q, r, s, t, p, K3());
  CGAL::internal::side_of_oriented_sphere(r, s, t, p, q, K3());
  CGAL::internal::side_of_oriented_sphere(s, t, p, q, r, K3());
  CGAL::internal::side_of_oriented_sphere(t, s, p, q, r, K3());
}

void
test_compare_squared_radius_3()
{
  // First test with random points.
  Point_3 p = my_rand_p3();
  Point_3 q = my_rand_p3();
  Point_3 r = my_rand_p3();
  Point_3 s = my_rand_p3();
  typedef std::pair<double,K4::FT> NT_pair;

  CGAL::internal::compare_squared_radius(p, q, r, s, NT_pair(0,0), K3());
  CGAL::internal::compare_squared_radius(p, q, r,    NT_pair(0,0), K3());
  CGAL::internal::compare_squared_radius(p, q,       NT_pair(0,0), K3());

  // Then with cospherical points (up to roundoff errors).
  p = sphere_rand_p3();
  q = sphere_rand_p3();
  r = sphere_rand_p3();
  s = sphere_rand_p3();

  CGAL::internal::compare_squared_radius(p, q, r, s, NT_pair(1,1), K3());
  //with three points
  double alpha_pqr = CGAL::to_double( K4().compute_squared_radius_3_object()(p.second,q.second,r.second) );
  CGAL::internal::compare_squared_radius(p, q, r, NT_pair(alpha_pqr,alpha_pqr), K3());
  //with 2 points
  double alpha_pq = CGAL::to_double( K4().compute_squared_radius_3_object()(p.second,q.second) );
  CGAL::internal::compare_squared_radius(p, q, NT_pair(alpha_pq,alpha_pq), K3());


  // Then with some perturbation.
  perturb(p, 1.0/(1<<25)/(1<<20)); // 2^-45

  CGAL::internal::compare_squared_radius(p, q, r, s, NT_pair(1,1), K3());
  CGAL::internal::compare_squared_radius(p, q, r,    NT_pair(alpha_pqr,alpha_pqr), K3());
  CGAL::internal::compare_squared_radius(p, q,       NT_pair(alpha_pq,alpha_pq), K3());
}


void compute_epsilons()
{
  K2::Orientation_2::compute_epsilon();
  K2::Orientation_3::compute_epsilon();
  K2::Side_of_oriented_circle_2::compute_epsilon();
  K2::Side_of_oriented_sphere_3::compute_epsilon();
}

int main(int argc, char **argv)
{
  int loops = (argc < 2) ? 2000 : std::atoi(argv[1]);
  int seed  = (argc < 3) ? CGAL::get_default_random().get_int(0, 1<<30)
                         : std::atoi(argv[2]);

  std::cout << "Initializing random generator with seed = " << seed
            << std::endl
            << "#loops = " << loops << " (can be changed on the command line)"
            << std::endl;

  CGAL::Random rnd(seed);
  r = &rnd;

  std::cout.precision(20);
  std::cerr.precision(20);

  #ifndef CGAL_NO_STATIC_FILTERS
  assert(  CGAL::Epick::Has_static_filters);
  #else
  assert( !CGAL::Epick::Has_static_filters);
  #endif

  #ifdef CGAL_DONT_USE_LAZY_KERNEL
    #ifdef CGAL_NO_STATIC_FILTERS
    assert(  !CGAL::Epeck::Has_static_filters);
    #else
    assert(  CGAL::Epeck::Has_static_filters);
    #endif
  #else
    #ifdef CGAL_LAZY_KERNEL_USE_STATIC_FILTERS_BY_DEFAULT
    assert(  CGAL::Epeck::Has_static_filters);
    #else
    assert( !CGAL::Epeck::Has_static_filters);
    #endif
  #endif

  assert(  K2::Has_static_filters);
  assert(! K4::Has_static_filters);

  std::cout << "ulp(1) = " << CGAL::internal::Static_filter_error::ulp() << std::endl;

  compute_epsilons();

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

  std::cout << "Testing Compare_squared_radius_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_compare_squared_radius_3();

  return 0;
}
