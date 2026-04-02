#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Double_grid_snap_rounding_traits_2.h>
#include <CGAL/vertical_slab_snap_rounding_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel                           Epeck;
typedef CGAL::Double_grid_snap_rounding_traits_2<Epeck, Epeck>                      Epeck_Traits;
typedef CGAL::Cartesian<CGAL::Exact_rational>                                       Rational_Kernel;
typedef CGAL::Double_grid_snap_rounding_traits_2<Rational_Kernel, Rational_Kernel>  Rational_Traits;

template <typename T>
struct minimal_output_iterator {
    // required nested type
    struct container_type {
        using value_type = T;
    };

    minimal_output_iterator& operator++() { return *this; }
    minimal_output_iterator operator++(int) { return *this; }
    minimal_output_iterator& operator*() { return *this; }
    minimal_output_iterator& operator=(const T&) { return *this; }
};

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
    CGAL::vertical_slab_snap_rounding_2(segs, std::back_inserter(out), CGAL::parameters::geom_traits(traits));
    assert(!CGAL::do_curves_intersect(out.begin(), out.end()));
    // Test minimal output iterator
    CGAL::vertical_slab_snap_rounding_2(segs, minimal_output_iterator< Segment_2 >(), CGAL::parameters::geom_traits(traits));
    CGAL::vertical_slab_snap_rounding_2(segs, minimal_output_iterator< std::vector<Point_2> >(), CGAL::parameters::geom_traits(traits));
  }
};

int main(/*int argc,char *argv[]*/)
{
  Test<Epeck_Traits>().fix_test();
  Test<Rational_Traits>().fix_test();
  return(0);
}
