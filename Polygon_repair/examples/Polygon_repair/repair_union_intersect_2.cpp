#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>

#include <CGAL/draw_multipolygon_with_holes_2.h>

#include <CGAL/Polygon_repair/repair.h>

#include <CGAL/IO/WKT.h>


#include <iostream>
#include <sstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point_2 = K::Point_2;
using Polygon_2 = CGAL::Polygon_2<K>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<K>;

int
main(int argc, char* argv[])
{
  Multipolygon_with_holes_2 pA;
  if (argc == 2) {
      {
          std::ifstream in(argv[1]);
          CGAL::IO::read_multi_polygon_WKT(in, pA);
      }
  } else {
          std::istringstream is("MULTIPOLYGON( ((0 0,  20 0, 20 20, 0 20), (1 1, 1 19, 19 19, 19 1) ) ,   (( 10 -2, 12 -2, 12 22, 10 22)) )");
          //std::istringstream is("MULTIPOLYGON( ((0 0, 2 0, 2 3, 0 3) ) )");    // (0.1 0.1, 0.1 0.4, 0.4 0.1)
          CGAL::IO::read_multi_polygon_WKT(is, pA);
  }

  Multipolygon_with_holes_2 mpwh = CGAL::Polygon_repair::repair(pA, CGAL::Polygon_repair::Union_rule());
  {
    std::ofstream out("union.wkt");
    CGAL::IO::write_multi_polygon_WKT(out, mpwh);
#ifdef CGAL_USE_BASIC_VIEWER
    CGAL::draw(mpwh);
#endif
  }
  mpwh = CGAL::Polygon_repair::repair(pA, CGAL::Polygon_repair::Intersection_rule());
  {
    std::ofstream out("intersection.wkt");
    CGAL::IO::write_multi_polygon_WKT(out, mpwh);

    CGAL::draw(mpwh);
  }

  {
    Polygon_2 pB;
    pB.push_back(Point_2(-1,-1));
    pB.push_back(Point_2(1,-1));
    pB.push_back(Point_2(1,1));
    pB.push_back(Point_2(-1,1));
    mpwh = CGAL::Polygon_repair::join(mpwh, pB);

    std::ofstream out("joinn.wkt");
    CGAL::IO::write_multi_polygon_WKT(out, mpwh);

    CGAL::draw(mpwh);
  }
  return 0;
}
