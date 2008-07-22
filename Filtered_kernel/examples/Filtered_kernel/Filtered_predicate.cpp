#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian_converter.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> FK;
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;
typedef CGAL::Cartesian_converter<K, EK> C2E;
typedef CGAL::Cartesian_converter<K, FK> C2F;

// Define my predicate, parameterized by a kernel.
template < typename K >
struct My_orientation_2
{
  typedef typename K::RT            RT;
  typedef typename K::Point_2       Point_2;

  typedef typename K::Orientation   result_type;

  result_type
  operator()(const Point_2 &p, const Point_2 &q, const Point_2 &r) const
  {
    RT prx = p.x() - r.x();
    RT pry = p.y() - r.y();
    RT qrx = q.x() - r.x();
    RT qry = q.y() - r.y();
    return CGAL::sign( prx*qry - qrx*pry );
  }
};

typedef CGAL::Filtered_predicate<My_orientation_2<EK>,
                                 My_orientation_2<FK>, C2E, C2F> Orientation_2;

int main()
{
  K::Point_2 p(1,2), q(2,3), r(3,4);
  Orientation_2 orientation;
  orientation(p, q, r);
  return 0;
}
