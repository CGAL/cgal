#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simplification_2/simplification_vw.h>
#include <CGAL/Simplification_2/simplification_ksbb.h>
#include <CGAL/Simplification_2/simplification_bmrs.h>
#include <CGAL/IO/WKT.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                       Polygon_with_holes_2;
typedef CGAL::Polyline_2<K>                                 Polyline_2;
typedef CGAL::StraightGeometryCollection<K>                 Collection;

void write_WKT(const std::string& filename, const Collection& collection)
{
  std::ofstream out(filename);
  out.precision(12);
  for(const auto& geometry : collection.geometries()) {
    std::visit(CGAL::overloaded{
        [&](const Point_2& p) { CGAL::IO::write_point_WKT(out, p); },
        [&](const Polygon_with_holes_2& p) { CGAL::IO::write_polygon_WKT(out, p); },
        [&](const Polyline_2& p) { CGAL::IO::write_linestring_WKT(out, p); },
        [](const auto&) {}},
      geometry);
  }
  std::cout << "  written to " << filename << "\n\n";
}

int main(int argc, char* argv[])
{
  // The default input is a small map: a west and an east region that share a border, a lake in the east region,
  // and a river that flows from the lake to the border.
  const std::string filename = (argc > 1) ? argv[1] : "data/map.wkt";
  // The fraction of the edges that remains after vw.
  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.3;

  std::ifstream in(filename);
  std::vector<Point_2> points;
  std::vector<Polyline_2> polylines;
  std::vector<Polygon_with_holes_2> polygons;
  if(!CGAL::IO::read_WKT(in, points, polylines, polygons)) {
    std::cerr << "Error: cannot read " << filename << std::endl;
    return 1;
  }

  // Geometries share a boundary or a vertex if they have the exact same points there.
  Collection input;
  for(const Polygon_with_holes_2& polygon : polygons)
      input.insert(polygon);
  for(const Polyline_2& polyline : polylines)
    input.insert(polyline);
  // Points are not simplified, but the simplified geometries keep them on the same side.
  for(const Point_2& point : points)
    input.insert(point);

  // Simplify with the VW vertex removal algorithm.
  CGAL::Simplification_VW vw(input);
  size_t inputComplexity = vw.complexity();
  vw.simplify(ratio * inputComplexity);
  Collection simplified_vw = vw.get_collection();
  write_WKT("simplified_vw.wkt", simplified_vw);

  // Simplify with the KSBB edge collapse algorithm.
  CGAL::Simplification_KSBB ksbb(input);
  ksbb.simplify(ratio * inputComplexity);
  Collection simplified_ksbb = ksbb.get_collection();
  write_WKT("simplified_ksbb.wkt", simplified_ksbb);

  // Make orthogonal schematization with the BMRS algorithm.
  CGAL::Simplification_BMRS bmrs(input, {.orientations=CGAL::regular_orientations<K>(2)});
  bmrs.simplify(ratio * inputComplexity); // note that here bmrs.complexity() is not the same as inputComplexity because of the orientation restriction step
  Collection simplified_bmrs = bmrs.get_collection();
  write_WKT("simplified_bmrs.wkt", simplified_bmrs);

  return 0;
}
