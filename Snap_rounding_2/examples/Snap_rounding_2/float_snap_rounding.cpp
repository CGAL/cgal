#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Float_snap_rounding_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Epick;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Curve_2                                       Segment_2;
typedef Kernel::Point_2                                         Point_2;
typedef Kernel::FT                                              FT;
typedef std::vector<Point_2 >                                   Polyline_2;

int main(int argc, char *argv[])
{
  std::vector< Segment_2 > segs;

  if(argc>1){
    std::cout << "Read segments in " << argv[1] << std::endl;
    std::ifstream in(argv[1]);

    int n;
    in >> n;
    segs.reserve(n);
    for (int i=0; i<n; ++i)
    {
      double x1,x2,y1,y2;
      in >> x1 >> y1 >> x2 >> y2;
      if(Point_2(x1,y1)!=Point_2(x2,y2))
        segs.emplace_back(Point_2(x1, y1), Point_2(x2, y2));
    }
  } else {
    //Example with a non-trivial rounding
    FT e(std::pow(2, -60));

    segs.emplace_back(Point_2(1-e, 1), Point_2(-1-e, -1+2*e));
    segs.emplace_back(Point_2(e/2, e/2), Point_2(1, -1));
    segs.emplace_back(Point_2(0, 2-e/2), Point_2(2, 0));
    segs.emplace_back(Point_2(0, 2-e/2), Point_2(-2+e, -4));
    segs.emplace_back(Point_2(-2, 2), Point_2(2, 2));
    segs.emplace_back(Point_2(7, 7), Point_2(7+e, 7+e));
    segs.emplace_back(Point_2(5, 7-e), Point_2(9, 7-e));
  }

  std::cout << "Computes the intersections and snaps the segments" << std::endl;
  std::vector< Segment_2> out;
  CGAL::Kernel_traits<std::remove_cv_t<std::iterator_traits<decltype(segs.begin())>::value_type>>::Kernel::toto toto;
  // CGAL::compute_snapped_subcurves_2<CGAL::Parallel_if_available_tag>(segs.begin(), segs.end(), out);
  // std::cout << "Does the output intersect: " << CGAL::do_curves_intersect(out.begin(), out.end()) << std::endl;
  std::cout << "Size of the output: " << out.size() << std::endl;

  std::string out_path=(argc>2)?argv[2]:"out.segs";
  std::cout << "Write the outputs in " << out_path << std::endl;
  std::ofstream outf(out_path);
  CGAL::Random r;

  outf << std::setprecision(17);
  outf << out.size() << std::endl;
  for (auto &seg: out)
  {
    double x1 = CGAL::to_double(seg.source().x());
    double y1 = CGAL::to_double(seg.source().y());
    double x2 = CGAL::to_double(seg.target().x());
    double y2 = CGAL::to_double(seg.target().y());
    outf << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
  }

  return 0;
}
