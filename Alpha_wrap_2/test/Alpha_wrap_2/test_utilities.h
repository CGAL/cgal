#ifndef CGAL_AW2_TEST_UTILITIES_2_H
#define CGAL_AW2_TEST_UTILITIES_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_repair/repair.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

// Reads polylines from .wkt or .obj
template <typename PolylineRange>
bool read_polylines(const std::string& filename,
                    PolylineRange& polylines)
{
  using Point_2 = typename std::decay<decltype(polylines.front().front())>::type;
  using Kernel = typename CGAL::Kernel_traits<Point_2>::type;

  using Points = std::vector<Point_2>;
  using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

  std::ifstream in(filename);
  if(!in)
    return false;

  const std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "wkt")
  {
    // read_WKT() reads ALL multi-linestrings whereas read_multilinestring() reads only the first one
    Points pts;
    Multipolygon mp;
    return CGAL::IO::read_WKT(in, pts, polylines, mp);
  }
  else if (ext == "obj")
  {
    std::vector<Point_2> points;
    std::vector<std::vector<std::size_t> > id_polylines;
    std::vector<std::vector<std::size_t> > unused_id_polygons;
    bool success = CGAL::IO::internal::read_OBJ(in, points, id_polylines, unused_id_polygons);
    if(!success)
      return false;

    for(const std::vector<std::size_t>& id_pl : id_polylines) {
      polylines.emplace_back();
      for(const std::size_t pid : id_pl) {
        polylines.back().push_back(points[pid]);
      }
    }

    return !polylines.empty();
  }
  return false;
}

template <typename Kernel>
void read_translate_and_write(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) {
    std::cerr << "Error: cannot open file '" << filename << "'" << std::endl;
    return;
  }

  using Point_2 = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;
  using Polyline = std::vector<Point_2>;
  using Polylines = std::vector<Polyline>;

  Polylines polylines;
  bool res = read_polylines(filename, polylines);
  if(!res) {
    std::cerr << "Error: cannot read polylines from file '" << filename << "'" << std::endl;
    return;
  }

  std::cout << "Read " << polylines.size() << " polylines" << std::endl;

  std::string out_filename = filename;
  std::size_t dot_pos = out_filename.find_last_of('.');
  if (dot_pos != std::string::npos) {
    out_filename = out_filename.substr(0, dot_pos) + "_translated.wkt";
  } else {
    out_filename += "_translated.wkt";
  }

  // compute the bbox of all polylines
  CGAL::Bbox_2 bbox = compute_bbox(polylines);

  // Translate all points by 20% of the vertical size downwards
  // and by 50% of the horizontal size to the left
  const double translate_x = -0.5 * (bbox.xmax() - bbox.xmin());
  const double translate_y = -0.2 * (bbox.ymax() - bbox.ymin());
  const Vector_2 translation(translate_x, translate_y);

  Polylines translated_polylines = polylines;
  for(Polyline& pl : translated_polylines) {
    for(Point_2& p : pl) {
      p = p + translation;
    }
  }

  polylines.insert(polylines.end(), translated_polylines.begin(), translated_polylines.end());

  std::ofstream out(out_filename);
  if (!out) {
    std::cerr << "Error: cannot write to file '" << out_filename << "'" << std::endl;
    return;
  }

  out.precision(std::numeric_limits<double>::max_digits10);
  CGAL::IO::write_multi_linestring_WKT(out, polylines);

  std::cout << "Wrote " << polylines.size() << " translated polylines to '" << out_filename << "'" << std::endl;
}

template <typename PolylineRange>
CGAL::Bbox_2 compute_bbox(const PolylineRange& pls)
{
  CGAL::Bbox_2 bbox;
  for (const auto& pl : pls) {
    for (const auto& p : pl) {
      bbox += p.bbox();
    }
  }
  return bbox;
}

template <typename MultipolygonWithHoles>
auto area(const MultipolygonWithHoles& mp)
{
  CGAL_precondition(CGAL::Polygon_repair::is_valid(mp));

  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using FT = typename Polygon_2::FT;

  FT total_area = FT(0);
  for(const Polygon_with_holes_2& pwh : mp.polygons_with_holes())
  {
    total_area += CGAL::abs(pwh.outer_boundary().area());
    for(const Polygon_2& hole : pwh.holes())
      total_area -= CGAL::abs(hole.area());
  }

  return total_area;
}

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_AW2_TEST_UTILITIES_2_H
