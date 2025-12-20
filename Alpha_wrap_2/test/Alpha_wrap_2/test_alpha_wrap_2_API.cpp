#define CGAL_AW2_TIMER

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Alpha_wrap_2/internal/validation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Random.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Alpha_wraps_2::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;
using Triangle_2 = Kernel::Triangle_2;

using Points = std::vector<Point_2>;

using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

// This is just to test the API, polygons are converted into polylines anyway

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  Multipolygon wrap;

  Polygon_2 poly;
  CGAL::alpha_wrap_2(poly, 1, wrap);
  CGAL::alpha_wrap_2(poly, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(poly, 1, 2, wrap);
  CGAL::alpha_wrap_2(poly, 1, 2, wrap, CGAL::parameters::all_default());

  Polygon_with_holes pwh;
  CGAL::alpha_wrap_2(pwh, 1, wrap);
  CGAL::alpha_wrap_2(pwh, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(pwh, 1, 2, wrap);
  CGAL::alpha_wrap_2(pwh, 1, 2, wrap, CGAL::parameters::all_default());

  Polygon_with_holes mp;
  CGAL::alpha_wrap_2(mp, 1, wrap);
  CGAL::alpha_wrap_2(mp, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(mp, 1, 2, wrap);
  CGAL::alpha_wrap_2(mp, 1, 2, wrap, CGAL::parameters::all_default());

  Points pts;
  CGAL::alpha_wrap_2(pts, 1, wrap);
  CGAL::alpha_wrap_2(pts, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(pts, 1, 2, wrap);
  CGAL::alpha_wrap_2(pts, 1, 2, wrap, CGAL::parameters::all_default());

  std::deque<Segment_2> sss;
  CGAL::alpha_wrap_2(sss, 1, wrap);
  CGAL::alpha_wrap_2(sss, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(sss, 1, 2, wrap);
  CGAL::alpha_wrap_2(sss, 1, 2, wrap, CGAL::parameters::all_default());

  std::vector<std::deque<Point_2> > mls;
  CGAL::alpha_wrap_2(mls, 1, wrap);
  CGAL::alpha_wrap_2(mls, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(mls, 1, 2, wrap);
  CGAL::alpha_wrap_2(mls, 1, 2, wrap, CGAL::parameters::all_default());

  std::list<Triangle_2> trs;
  CGAL::alpha_wrap_2(trs, 1, wrap);
  CGAL::alpha_wrap_2(trs, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(trs, 1, 2, wrap);
  CGAL::alpha_wrap_2(trs, 1, 2, wrap, CGAL::parameters::all_default());

  std::vector<Point_2> pos;
  std::vector<std::array<std::size_t, 3> > faces;
  CGAL::alpha_wrap_2(pos, faces, 1, wrap);
  CGAL::alpha_wrap_2(pos, faces, 1, wrap, CGAL::parameters::all_default());
  CGAL::alpha_wrap_2(pos, faces, 1, 2, wrap);
  CGAL::alpha_wrap_2(pos, faces, 1, 2, wrap, CGAL::parameters::all_default());

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
