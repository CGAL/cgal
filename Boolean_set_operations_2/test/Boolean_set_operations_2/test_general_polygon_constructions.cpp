#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2/Polygon_conversions.h>

template <typename ArrTraits, typename Kernel, typename Container>
CGAL::Polygon_2<Kernel, Container>
my_convert_polygon_back(CGAL::General_polygon_2<ArrTraits>& gpgn,
                        CGAL::Polygon_2<Kernel, Container>&)
{ return CGAL::convert_polygon_back<Kernel, Container>(gpgn); }

template <typename ArrTraits, typename Kernel, typename Container>
CGAL::Polygon_with_holes_2<Kernel, Container>
my_convert_polygon_back(CGAL::General_polygon_with_holes_2
                        <CGAL::General_polygon_2<ArrTraits>>& gpwh,
                        CGAL::Polygon_with_holes_2<Kernel, Container>&)
{ return CGAL::convert_polygon_back<Kernel, Container>(gpwh); }

int main() {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Point_2 = Kernel::Point_2;
  using Segment_2 = Kernel::Segment_2;
  using Segment_traits = CGAL::Arr_segment_traits_2<Kernel>;
  using Polyline_traits = CGAL::Arr_polyline_traits_2<Segment_traits>;
  using Gps_traits = CGAL::Gps_traits_2<Polyline_traits>;

  using General_pgn = Gps_traits::General_polygon_2;
  using General_pwh = Gps_traits::General_polygon_with_holes_2;

  using X_monotone_curve_2 = Polyline_traits::X_monotone_curve_2;

  using Pgn = CGAL::Polygon_2<Kernel>;
  using Pwh = CGAL::Polygon_with_holes_2<Kernel>;

  std::vector<Point_2> points = {
    Point_2(0, 0),
    Point_2(1, 0),
    Point_2(2, 1),
    Point_2(1, 1),
    Point_2(0, 0.5)
  };

  auto point_to_segment = [&](const Point_2& p) -> Segment_2 {
    const Point_2* it = &p;
    if (it == &points.back()) return Segment_2(p, points[0]);
    return Segment_2(p, *(++it));
  };

  Gps_traits gtraits;
  const Polyline_traits& ptraits(gtraits);
  auto ctr = ptraits.construct_curve_2_object();
  auto make_x_mtn = ptraits.make_x_monotone_2_object();
  auto eql = gtraits.equal_2_object();

  // Case 1: Segment-based GPS from Segment_2 range with transform iterator
  CGAL::General_polygon_2<Segment_traits> gpgn1
    (boost::make_transform_iterator(points.begin(), point_to_segment),
     boost::make_transform_iterator(points.end(), point_to_segment));
  std::cout << "gpgn1: " << gpgn1 << std::endl;

  // Case 2: Polyline-based GPS from polyline of points
  auto curve2 =
    ctr(boost::range::join(CGAL::make_range(points.begin(), points.end()),
                           CGAL::make_single(*points.begin())));
  General_pgn gpgn2;
  using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;
  make_x_mtn(curve2,
             boost::make_function_output_iterator
             ([&](const Make_x_monotone_result& obj)
              { gpgn2.push_back(*(std::get_if<X_monotone_curve_2>(&obj))); }));
  std::cout << "gpgn2: " << gpgn2 << std::endl;

  // Case 3: Polyline-based GPS from polyline of segments
  auto curve3 =
    ctr(boost::make_transform_iterator(points.begin(), point_to_segment),
        boost::make_transform_iterator(points.end(), point_to_segment));
  General_pgn gpgn3;
  make_x_mtn(curve3,
             boost::make_function_output_iterator
             ([&](const Make_x_monotone_result& obj)
              { gpgn3.push_back(*(std::get_if<X_monotone_curve_2>(&obj))); }));
  std::cout << "gpgn3: " << gpgn3 << std::endl;

  if (! eql(gpgn2, gpgn3)) {
    std::cerr << "Construction 1 failed\n";
    return EXIT_FAILURE;
  }

  // Case 4: Polygon from From points
  Pgn pgn1(points.begin(), points.end());
  std::cout << "pgn1: " << pgn1 << std::endl;

  // Case 5: Polyline-based GPS from conversion
  General_pgn gpgn4 = CGAL::convert_polygon(pgn1, ptraits);
  std::cout << "gpgn4: " << gpgn4 << std::endl;

  if (! eql(gpgn4, gpgn3)) {
    std::cerr << "Construction 2 failed\n";
    return EXIT_FAILURE;
  }

  // Case 5: Polyline-based GPS from conversion
  Pgn pgn2 = my_convert_polygon_back(gpgn4, pgn1);
  std::cout << "pgn2: " << pgn2 << std::endl;
  if (pgn1 != pgn2) {
    std::cerr << "Construction 2 failed\n";
    return EXIT_FAILURE;
  }

  Pwh pwh1(pgn1);
  std::cout << "pwh1: " << pwh1 << std::endl;

  General_pwh gpwh1(gpgn4);
  std::cout << "gpwh1: " << gpwh1 << std::endl;

  General_pwh gpwh2 = convert_polygon(pwh1, ptraits);
  std::cout << "gpwh2: " << gpwh2 << std::endl;

  if (! eql(gpwh1, gpwh2)) {
    std::cerr << "Construction 3 failed\n";
    return EXIT_FAILURE;
  }

  Pwh pwh2 = my_convert_polygon_back(gpwh1, pwh1);
  std::cout << "pwh2: " << pwh2 << std::endl;
  if (pwh1 != pwh2) {
    std::cerr << "Construction 4 failed\n";
    return EXIT_FAILURE;
  }

  Pwh pwh3;
  CGAL::Oneset_iterator<Pwh> iterator(pwh3);
  auto converter = convert_polygon_back(iterator, pwh1);
  *converter++ = gpwh1;
  std::cout << "pwh3: " << pwh3 << std::endl;

  if (pwh1 != pwh3) {
    std::cerr << "Construction 5 failed\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
