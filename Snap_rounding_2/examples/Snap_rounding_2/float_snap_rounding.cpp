#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Float_snap_rounding_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Curve_2                                       Segment_2;
typedef Kernel::Point_2                          Point_2;
typedef Kernel::FT                               FT;

int main()
{
  std::vector<Point_2> pts;
  std::vector< Segment_2 > segs;

  Kernel::FT e(std::pow(2, -60));

  segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
  segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
  segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
  segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
  segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
  segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));

  std::cout << "Input" << std::endl;
  for(auto &seg: segs){
    std::cout << seg.source() << " -> " << seg.target() << std::endl;
  }
  std::cout << "\n\n" << std::endl;

  std::vector< std::vector<Point_2 > > out;
  CGAL::double_snap_rounding_2(segs.begin(), segs.end(), out);

  std::cout << "Output" << std::endl;
  for(auto &poly: out){
    std::cout << poly[0];
    for(size_t i=1; i!=poly.size(); ++i)
        std::cout << " -> " << poly[i];
    std::cout << std::endl;
  }

  return 0;
}
