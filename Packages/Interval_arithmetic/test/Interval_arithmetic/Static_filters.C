
#define CGAL_PROFILE

#include <CGAL/Random.h>
#include <CGAL/Static_filters.h>

typedef CGAL::Simple_cartesian<double> K;

CGAL::Random R;

K::Point_3 prand()
{
  return K::Point_3(R.get_double(), R.get_double(), R.get_double());
}

K::Point_2 prand2()
{
  return K::Point_2(R.get_double(), R.get_double());
}

CGAL::Orientation ooo;
CGAL::Oriented_side os;

int main()
{
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

  for (int i=0; i<100; ++i) {
    K::Point_3 p(prand()), q(prand()), r(prand()), s(prand()), t(prand());
    K::Point_2 p2(prand2()), q2(prand2()), r2(prand2()), s2(prand2());
    k.register_object(p);
    k.register_object(q);
    k.register_object(r);
    k.register_object(s);
    k.register_object(t);
    k.register_object(p2);
    k.register_object(q2);
    k.register_object(r2);
    k.register_object(s2);
    // CGAL::Protect_FPU_rounding<false> Z;
    for (int j=0; j<1000; ++j) {
      // assert( my_o(p, q, r, s) == ore(pe, qe, re, se) );
      ooo = my_o3(p, q, r, s);      //  2.62 s -> 2.22 s
      ooo = my_o2(p2, q2, r2);      //  2.62 s -> 2.22 s
      os = my_c2(p2, q2, r2, s2);
      os = my_s3(p, q, r, s, t);
      ooo = my_co3(p, q, r);
      ooo = my_co3(p, q, r, s);

      // ooo = ore(pe, qe, re, se);   // 15.07 s  (!prot:13.5) static: 4.35 s (!p:3 -> 1.45)
      // ooo = orf(pf, qf, rf, sf);   // 0.68 s -> 0.85 s
      // ooo = inexact_o(p, q, r, s); //  1.66 s -> 0.20 s
    }
  }

  std::cerr << "Number of IA failures = "
            << CGAL::Interval_base::number_of_failures << std::endl;

  return 0;
}
