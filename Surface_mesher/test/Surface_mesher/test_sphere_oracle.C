#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesher/Sphere_oracle_3.h>

#include <CGAL/point_generators_3.h>

#include <functional>
#include <iostream>

#include <stdlib.h>

#include <boost/array.hpp>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef K::Point_3 Point_3;

typedef K::Segment_3 Segment_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;

typedef K::Sphere_3 Sphere_3;

typedef K::FT FT;

class Sphere {
private:
  const FT r; 
  const Point_3 c;
public:
  Sphere(FT rayon, Point_3 c = CGAL::ORIGIN) : r(rayon), c(c) {}

  FT operator()(const Point_3& p) const
  {
    const FT x = p.x()-c.x();
    const FT y = p.y()-c.y();
    const FT z = p.z()-c.z();
    return x*x+y*y+z*z-r*r;
  }
};

typedef CGAL::Implicit_surface_3<K, Sphere> Implicit_sphere;

typedef CGAL::Surface_mesher::Implicit_surface_oracle_3<K, Implicit_sphere> Implicit_oracle;
typedef CGAL::Surface_mesher::Sphere_oracle_3<K> Sphere_oracle;

typedef Implicit_oracle::Intersect_3 Implicit_intersect_3;
typedef Sphere_oracle::Intersect_3 Sphere_intersect_3;

struct Too_close_to_sphere : public std::unary_function<Point_3, bool>
{
  const K::FT radius;
  const K::FT precision;
  
  Too_close_to_sphere(const K::FT& radius, K::FT precision = 0.1)
    : radius(radius), precision(precision)
  {
  }
  
  template <typename Point_iterator>
  bool operator()(Point_iterator it) const
  {
    const Point_3& p = *it;
    const K::FT sq_distance = CGAL::squared_distance(Point_3(0, 0, 0), p);
    const K::FT sq_radius = radius*radius;

    return sq_radius * (1.-precision) <= sq_distance && 
      sq_distance <= (1.+precision) * sq_radius;
  }
};

struct In_sphere : public std::unary_function<Point_3, bool>
{
  const K::FT radius;
  const K::FT precision;
  
  In_sphere(const K::FT& radius, K::FT precision = 0.1)
    : radius(radius), precision(precision)
  {
  }
  
  template <typename Point_iterator>
  bool operator()(Point_iterator it) const
  {
    const Point_3& p = *it;
    const K::FT sq_distance = CGAL::squared_distance(Point_3(0, 0, 0), p);
    const K::FT sq_radius = radius*radius;

    return sq_distance < sq_radius * (1.+precision);

  }
};

struct Out_of_sphere : public std::unary_function<Point_3, bool>
{
  const K::FT radius;
  const K::FT precision;
  
  Out_of_sphere(const K::FT& radius, K::FT precision = 0.1)
    : radius(radius), precision(precision)
  {
  }
  
  template <typename Point_iterator>
  bool operator()(Point_iterator it) const
  {
    const Point_3& p = *it;
    const K::FT sq_distance = CGAL::squared_distance(Point_3(0, 0, 0), p);
    const K::FT sq_radius = radius*radius;

    return sq_distance > (1.-precision) * sq_radius;
  }
};
const K::FT radius = 2;

int main(int, char**)
{
  Implicit_intersect_3 implicit_intersect = Implicit_oracle().intersect_3_object();
  Sphere_intersect_3 intersect = Sphere_oracle().intersect_3_object();

  int exit = EXIT_SUCCESS;


  // first test
  std::cerr << "First test...\n";

  Sphere_3 test_1(Point_3(2,3,4), 5*5);

  std::cerr << "sphere: " << test_1 << "\n";

  boost::array<Segment_3, 5> segs;
  segs[0] = Segment_3(Point_3(-4, 3, 4), Point_3(8, 3, 4));
  segs[1] = Segment_3(Point_3(2, 9, 4), Point_3(2, -3, 4));
  segs[2] = Segment_3(Point_3(3, 3, 4), Point_3(8, 10, -4));
  segs[3] = Segment_3(Point_3(4, 3, 4), Point_3(8, -10, 4));
  segs[4] = Segment_3(Point_3(2, 3, -2), Point_3(2, 3, 10));

  for(unsigned int i = 0; i < segs.size(); ++i)
  {
    std::cerr << "segment " << i << ": " << segs[i] << "\n";

    CGAL::Object o = intersect(test_1, segs[i]);
    if( const Point_3* p_pt = CGAL::object_cast<Point_3>(&o) )
    {
      std::cerr << "intersection: " << *p_pt 
                << " (distance=" << CGAL::sqrt(CGAL::squared_distance(test_1.center(), *p_pt))
                << ")\n";
    }
    else
      exit = EXIT_FAILURE;

    Point_3 a = segs[i].source();
    Point_3 b = segs[i].target();
    intersect.clip_segment(test_1, a, b);
    std::cerr << "clipped segment: " << a << ", " << b << "\n"
                << " (distance: " << CGAL::sqrt(CGAL::squared_distance(test_1.center(), a))
                << ", " << CGAL::sqrt(CGAL::squared_distance(test_1.center(), b))
                << ")\n\n";
  }

  // second test
  std::cerr << "Second test...\n";

  typedef CGAL::Random_points_in_cube_3<Point_3> Random_points_in_cube;

  Random_points_in_cube in_cube(10.);

  CGAL::Filter_iterator<Random_points_in_cube, 
    Too_close_to_sphere > points(Random_points_in_cube(1),
                                 Too_close_to_sphere(radius),
                                 Random_points_in_cube(radius*1.1));

  CGAL::Filter_iterator<Random_points_in_cube, 
    Out_of_sphere > points_out(Random_points_in_cube(1),
                               Out_of_sphere(radius),
                               Random_points_in_cube(radius*1.1));

  CGAL::Filter_iterator<Random_points_in_cube, 
    In_sphere > points_in(Random_points_in_cube(1),
                          In_sphere(radius),
                          Random_points_in_cube(radius*1.1));

  Point_3 center = *in_cube++;

  const Sphere_3 bounding_sphere(center, radius*radius*4);

  Implicit_sphere implicit_sphere(Sphere(radius, // not squared radius,
                                                 // see class ::Sphere
                                         center),
                                  bounding_sphere);

  Sphere_3 sphere(center, radius*radius);
  
  unsigned int number_of_segment_intersections = 0;
  unsigned int number_of_ray_intersections = 0;
  unsigned int number_of_line_intersections = 0;
  unsigned int number_of_impl_segment_intersections = 0;
  unsigned int number_of_impl_ray_intersections = 0;
  unsigned int number_of_impl_line_intersections = 0;
  unsigned int number_of_problems = 0;
  unsigned int number_of_both_outside = 0;

  std::cerr << "sphere: " << sphere << std::endl;
  std::cerr << "bounding sphere: " << bounding_sphere << std::endl;

  for(int n = 0; n < 1500; ++n)
  {
    const Point_3 a = *points++ + ( center - CGAL::ORIGIN );
    const Point_3 b = *points++ + ( center - CGAL::ORIGIN );

    const Segment_3 s = Segment_3(a, b);
    const Ray_3 r = Ray_3(a, b);
    const Line_3 l = Line_3(a, b);

    CGAL::Object inter_segment_impl = implicit_intersect(implicit_sphere, s);
    CGAL::Object inter_ray_impl =     implicit_intersect(implicit_sphere, r);
    CGAL::Object inter_line_impl =    implicit_intersect(implicit_sphere, l);
    CGAL::Object inter_segment =      intersect(sphere, s);
    CGAL::Object inter_ray =          intersect(sphere, r);
    CGAL::Object inter_line =         intersect(sphere, l);

    const bool test_segment_impl = ! inter_segment_impl.is_empty();
    const bool test_ray_impl =     ! inter_ray_impl .is_empty();
    const bool test_line_impl =    ! inter_line_impl.is_empty();
    const bool test_segment =      ! inter_segment.is_empty();
    const bool test_ray =          ! inter_ray.is_empty();
    const bool test_line =         ! inter_line.is_empty();

    const bool a_in_sphere = sphere.has_on_bounded_side(a);
    const bool b_in_sphere = sphere.has_on_bounded_side(b);

    const bool both_outside = !(a_in_sphere || b_in_sphere);

    if(both_outside) 
      ++number_of_both_outside;

    bool problem = false;

    if( test_segment_impl )
    {
      ++number_of_impl_segment_intersections;
    }

    if( test_ray_impl )
    {
      ++number_of_impl_ray_intersections;
    }

    if( test_line_impl )
    {
      ++number_of_impl_line_intersections;
    }

    if( test_segment )
    {
      ++number_of_segment_intersections;
    }

    if( test_ray )
    {
      ++number_of_ray_intersections;
    }

    if( test_line )
    {
      ++number_of_line_intersections;
    }
    if( a_in_sphere != b_in_sphere )
    {
      if( ! test_segment )
      {
        std::cerr << "Should be a segment intersection.\n";
        problem = true;
      }
      if( ! test_ray )
      {
        std::cerr << "Should be a ray intersection.\n";
        problem = true;
      }
      if( ! test_line )
      {
        std::cerr << "Should be a line intersection.\n";
        problem = true;
      }
    }

    if( !both_outside && test_segment_impl != test_segment)
    {
      std::cerr << "Problem with segment intersection.\n";
      if(test_segment_impl)
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_segment_impl);
        std::cerr << "inter_segment_impl=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      else
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_segment);
        std::cerr << "inter_segment=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      problem = true;
    }

    if( !both_outside &&  test_ray_impl != test_ray)
    {
      std::cerr << "Problem with ray intersection.\n";
      if(test_ray_impl)
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_ray_impl);
        std::cerr << "inter_ray_impl=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      else
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_ray);
        std::cerr << "inter_ray=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      problem = true;
    }

    if( !both_outside &&  test_line_impl != test_line)
    {
      std::cerr << "Problem with line intersection.\n";
      if(test_line_impl)
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_line_impl);
        std::cerr << "inter_line_impl=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      else
      {
        const Point_3* p;
        p = CGAL::object_cast<Point_3>(&inter_line);
        std::cerr << "inter_line=" << *p
                  << "\ndistance: " << CGAL::sqrt(CGAL::squared_distance(center, *p))
                  << "\n";
      }
      problem = true;
    }

    if( test_segment )
      if( ! test_ray ) 
      {
        std::cerr << "Segment intersection, but no ray intersection! \n";
        problem = true;
      }

    if( test_segment || test_ray )
      if( ! test_line ) 
      {
        std::cerr << "Segment or ray intersection, but no line intersection! \n";
        problem = true;
      }
    
    if(problem)
    {
      ++number_of_problems;
      std::cerr << "segment: " << s
                << "\n(Distances: " << CGAL::sqrt(CGAL::squared_distance(center, a)) << " - "
                << CGAL::sqrt(CGAL::squared_distance(center, b)) << ")\n";

      Point_3 c = a;
      Point_3 d = b;
      intersect.clip_line(sphere, Line_3(a, d), c, d);
      std::cerr << "clipped line: " << Segment_3(c, d)
                << "\n(Distances: " << CGAL::sqrt(CGAL::squared_distance(center, c)) << " - "
                << CGAL::sqrt(CGAL::squared_distance(center, d)) << ")\n";

      intersect.clip_ray(sphere, Ray_3(a, d), c, d);
      std::cerr << "clipped ray: " << Segment_3(c, d)
                << "\n(Distances: " << CGAL::sqrt(CGAL::squared_distance(center, c)) << " - "
                << CGAL::sqrt(CGAL::squared_distance(center, d)) << ")\n";

      c = a;
      d = b;
      intersect.clip_segment(sphere, c, d);
      std::cerr << "clipped segment: " << Segment_3(c, d)
                << "\n(Distances: " << CGAL::sqrt(CGAL::squared_distance(center, c)) << " - "
                << CGAL::sqrt(CGAL::squared_distance(center, d)) << ")\n\n";
      exit = EXIT_FAILURE;
    }
  }

  std::cerr
    <<   "Number of segment intersections: " << number_of_segment_intersections
    << "\nNumber of ray intersections: " << number_of_ray_intersections
    << "\nNumber of line intersections: " << number_of_line_intersections
    << "\nNumber of impl_segment intersections: " << number_of_impl_segment_intersections
    << "\nNumber of impl_ray intersections: " << number_of_impl_ray_intersections
    << "\nNumber of impl_line intersections: " << number_of_impl_line_intersections
    << "\nNumber of both outside: " << number_of_both_outside
    << "\nNumber of problems: " << number_of_problems
    << "\n";
  return exit;
}
