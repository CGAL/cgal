
#define CGAL_PROFILE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Static_filters.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Random_points_in_square_2<K::Point_2> Rand_2;

CGAL::Orientation   oo;
CGAL::Oriented_side os;
CGAL::Bounded_side  bs;

int main()
{
  CGAL::Random_points_in_square_2<K::Point_2> Rand_2(1000.0);
  CGAL::Random_points_in_cube_3<K::Point_3>   Rand_3(1000.0);
  CGAL::Static_filters<K> k;

  const CGAL::Static_filters<K>::Orientation_3 & my_o3 =
        k.orientation_3_object();
  const CGAL::Static_filters<K>::Orientation_2 & my_o2 =
        k.orientation_2_object();
  const CGAL::Static_filters<K>::Side_of_oriented_circle_2 & my_c2 =
        k.side_of_oriented_circle_2_object();
  const CGAL::Static_filters<K>::Side_of_oriented_sphere_3 & my_s3 =
        k.side_of_oriented_sphere_3_object();
  const CGAL::Static_filters<K>::Coplanar_orientation_3 & my_co3 =
        k.coplanar_orientation_3_object();
  const CGAL::Static_filters<K>::Coplanar_side_of_bounded_circle_3 & my_cobc3 =
        k.coplanar_side_of_bounded_circle_3_object();

  // my_cobc3.cir_3();

  for (int i=0; i<100; ++i) {
    K::Point_3 p = *Rand_3++; k.register_object(p);
    K::Point_3 q = *Rand_3++; k.register_object(q);
    K::Point_3 r = *Rand_3++; k.register_object(r);
    K::Point_3 s = *Rand_3++; k.register_object(s);
    K::Point_3 t = *Rand_3++; k.register_object(t);

    K::Point_2 p2 = *Rand_2++; k.register_object(p2);
    K::Point_2 q2 = *Rand_2++; k.register_object(q2);
    K::Point_2 r2 = *Rand_2++; k.register_object(r2);
    K::Point_2 s2 = *Rand_2++; k.register_object(s2);

    for (int j=0; j<1000; ++j) {
      oo = my_o3(p, q, r, s);
      oo = my_o2(p2, q2, r2);
      os = my_c2(p2, q2, r2, s2);
      os = my_s3(p, q, r, s, t);
      oo = my_co3(p, q, r);
      oo = my_co3(p, q, r, s);
      bs = my_cobc3(p, q, r, s);
    }
  }

  std::cerr << "Number of IA failures = "
            << CGAL::Interval_base::number_of_failures << std::endl;

  return 0;
}
