
#define CGAL_PROFILE

#include <CGAL/Random.h>
#include <CGAL/Fixed_precision_nt.h>
#include <CGAL/Static_filters.h>

typedef CGAL::Simple_cartesian<CGAL::Fixed_precision_nt> Kf;
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Simple_cartesian<CGAL::Filtered_exact<double, CGAL::MP_Float> > //, CGAL::Dynamic, false> >
// typedef CGAL::Simple_cartesian<CGAL::Filtered_exact<double, CGAL::MP_Float, CGAL::Static, false> >
        KE;

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

int main()
{
  CGAL::Static_filters<K> k;
  const CGAL::Static_filters<K>::Orientation_3 & my_o=k.orientation_3_object();
  const CGAL::Static_filters<K>::Orientation_2 & my_o2=k.orientation_2_object();
  // K::Orientation_3 inexact_o;
  // KE::Orientation_3 ore;
  // Kf::Orientation_3 orf;

#ifdef CGAL_IA_NEW_FILTERS
  CGAL::Static_Filtered_orientationC3_12::new_bound(1.0);
#endif
  CGAL::Fixed_precision_nt::init(1.0);
  std::cerr << "*** bound used by Fixed : " << CGAL::Fixed_Or1 << std::endl;
  // This bound is slightly better than epsilon_2.  It's OK, since for the
  // Fixed, the initial subtraction fits exactly in a double.
  // Maybe we can slightly improve the bound, but not much, it seems.

  for (int i=0; i<100; ++i) {
    K::Point_3 p(prand()), q(prand()), r(prand()), s(prand());
    K::Point_2 p2(prand2()), q2(prand2()), r2(prand2());
    // KE::Point_3 pe(p.x(), p.y(), p.z()), qe(q.x(), q.y(), q.z()),
    //             re(r.x(), r.y(), r.z()), se(s.x(), s.y(), s.z());
    // Kf::Point_3 pf(p.x(), p.y(), p.z()), qf(q.x(), q.y(), q.z()),
    //             rf(r.x(), r.y(), r.z()), sf(s.x(), s.y(), s.z());
    // CGAL::Protect_FPU_rounding<false> Z;
    for (int j=0; j<100000; ++j) {
      // assert( my_o(p, q, r, s) == ore(pe, qe, re, se) );
      ooo = my_o(p, q, r, s);      //  2.62 s -> 2.22 s
      ooo = my_o2(p2, q2, r2);      //  2.62 s -> 2.22 s
      // ooo = ore(pe, qe, re, se);   // 15.07 s  (!prot:13.5) static: 4.35 s (!p:3 -> 1.45)
      // ooo = orf(pf, qf, rf, sf);   // 0.68 s -> 0.85 s
      // ooo = inexact_o(p, q, r, s); //  1.66 s -> 0.20 s
    }
  }

  std::cerr << "Number of IA failures = " << CGAL::Interval_base::number_of_failures << std::endl;

  return 0;
}
