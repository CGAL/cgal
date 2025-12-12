// #define DOUBLE_2D_SNAP_VERBOSE
// #define DOUBLE_2D_SNAP_FULL_VERBOSE
#define BENCH_AND_VERBOSE_FLOAT_SNAP_ROUNDING_2

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Float_snap_rounding_traits_2.h>
#include <CGAL/Float_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel                      Epeck;
typedef CGAL::Cartesian<CGAL::Exact_rational>                                 Rational_Kernel;
typedef CGAL::Float_snap_rounding_traits_2<Rational_Kernel, Rational_Kernel>  Rational_Traits;

template<class Traits>
struct Test{
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Point_2   Point_2;

  void fix_test(){
    std::vector<Segment_2> segs;
    segs.emplace_back(Point_2(1, 1), Point_2(-1, -1));
    segs.emplace_back(Point_2(0, 0), Point_2(1, -1));
    segs.emplace_back(Point_2(0, 2), Point_2(2, 0));
    segs.emplace_back(Point_2(0, 2), Point_2(-2, -4));
    segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
    segs.emplace_back(Point_2(5, 7), Point_2(9, 7));
    segs.emplace_back(Point_2(1,1), Point_2(3,1));
    segs.emplace_back(Point_2(1,2), Point_2(3,0));
    std::vector<Segment_2> out;
    Traits traits;
    CGAL::compute_snapped_subcurves_2(segs.begin(), segs.end(), std::back_inserter(out), CGAL::parameters::geom_traits(traits));
    assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
  }
};

int main(/*int argc,char *argv[]*/)
{
  Test<Rational_Traits>().fix_test();
  return(0);
}
