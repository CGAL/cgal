#ifndef CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H
#define CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H

// STL includes.
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// CGAL includes.
#include <CGAL/IO/io.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h>

template<typename GeomTraits>
struct Saver {

public:
  using Traits = GeomTraits;
  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Polyline = std::vector<Point_3>;

  Saver() {
    out.precision(20);
  }

  inline std::string data() const {
    return out.str();
  }

  void export_segments(
    const std::vector<Segment_2>& segments,
    const std::string path,
    const FT) {

    std::vector<Polyline> polylines(segments.size());
    for (std::size_t i = 0; i < segments.size(); ++i) {
      const auto& s = segments[i].source();
      const auto& t = segments[i].target();

      polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
      polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
    }
    export_polylines(polylines, path);
  }

  void export_group(
    const std::vector<Segment_2>& segments,
    const std::vector<std::size_t>& group,
    const std::string path,
    const FT) {

    const FT stub = FT(0);
    std::vector<Segment_2> edges;
    for (const std::size_t seg_index : group) {
      edges.push_back(segments[seg_index]);
    }
    export_segments(edges, path, stub);
  }

  void export_closed_contour(
    const std::vector<Point_2>& contour,
    const std::string path,
    const FT) {

    if (contour.size() == 0) {
      return;
    }

    const FT stub = FT(0);
    std::vector<Segment_2> segments;
    const std::size_t n = contour.size();
    segments.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& p = contour[i];
      const auto& q = contour[ip];
      segments.push_back(Segment_2(p, q));
    }
    export_segments(segments, path, stub);
  }

  void export_open_contour(
    const std::vector<Point_2>& contour,
    const std::string path,
    const FT) {

    if (contour.size() == 0) {
      return;
    }

    const FT stub = FT(0);
    std::vector<Segment_2> segments;
    const std::size_t n = contour.size();
    segments.reserve(n - 1);

    for (std::size_t i = 0; i < n - 1; ++i) {
      const std::size_t ip = i + 1;

      const auto& p = contour[i];
      const auto& q = contour[ip];
      segments.push_back(Segment_2(p, q));
    }
    export_segments(segments, path, stub);
  }

  void export_eps_segments(
    const std::vector<Segment_2>& input,
    const std::string path,
    FT scale) {

    if (input.size() == 0) return;
    clear();

    // Compute barycenter.
    Point_2 b;
    compute_barycenter(input, b);

    // Translate segments.
    std::vector<Segment_2> segments;
    translate_segments(input, b, segments);

    // Compute bounding box.
    Point_2 minb, maxb;
    compute_bounding_box(segments, minb, maxb);

    // Estimate eps parameters.
    const FT length = static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(minb, maxb))));
    if (length < FT(10) && scale == FT(1)) scale *= FT(1000);

    // const FT radius = FT(1);
    const FT line_width = FT(1);
    const bool dashed = false;

    // Set eps header.
    set_eps_header(
      minb.x() * scale, minb.y() * scale,
      maxb.x() * scale, maxb.y() * scale,
      "segments");

    // Start private namespace.
    out << "0 dict begin gsave" << std::endl << std::endl;

    // Draw segments.
    for (const auto& segment : segments) {
      add_eps_segment(segment, scale, line_width, dashed);
      // add_eps_disc(segment.source(), radius, scale);
      // add_eps_disc(segment.target(), radius, scale);
    }

    // Finish private namespace.
    out << "grestore end" << std::endl << std::endl;
    out << "%%EOF" << std::endl;
    save(path + ".eps");
  }

  void export_eps_group(
    const std::vector<Segment_2>& segments,
    const std::vector<std::size_t>& group,
    const std::string path) {

    const FT stub = FT(0);
    std::vector<Segment_2> edges;
    for (const std::size_t seg_index : group) {
      edges.push_back(segments[seg_index]);
    }
    export_eps_segments(edges, path, stub);
  }

  void export_eps_closed_contour(
    const std::vector<Point_2>& contour,
    const std::string path,
    FT scale) {

    if (contour.size() == 0) {
      return;
    }

    std::vector<Segment_2> segments;
    const std::size_t n = contour.size();
    segments.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& p = contour[i];
      const auto& q = contour[ip];
      segments.push_back(Segment_2(p, q));
    }
    export_eps_segments(segments, path, scale);
  }

  void export_eps_open_contour(
    const std::vector<Point_2>& contour,
    const std::string path,
    FT scale) {

    if (contour.size() == 0) {
      return;
    }

    std::vector<Segment_2> segments;
    const std::size_t n = contour.size();
    segments.reserve(n - 1);

    for (std::size_t i = 0; i < n - 1; ++i) {
      const std::size_t ip = i + 1;

      const auto& p = contour[i];
      const auto& q = contour[ip];
      segments.push_back(Segment_2(p, q));
    }
    export_eps_segments(segments, path, scale);
  }

private:
  std::stringstream out;

  void clear() {
    out.str(std::string());
  }

  void save(
    const std::string path) const {

    std::ofstream file(path.c_str(), std::ios_base::out);
    CGAL::IO::set_ascii_mode(file);
    if (!file) {
      std::cerr <<
        "Error: cannot save the file: " << path << std::endl; return;
    }

    file << data() << std::endl; file.close();
    std::cout <<
      "* data are saved in " << path << std::endl;
  }

  void export_polylines(
    const std::vector<Polyline>& polylines,
    const std::string path) {

    if (polylines.size() == 0) {
      return;
    }

    clear();
    for (std::size_t i = 0; i < polylines.size(); ++i) {
      const auto& polyline = polylines[i];

      out << polyline.size() << " ";
      for (std::size_t j = 0; j < polyline.size(); ++j)
        out << polyline[j] << " ";
      out << std::endl;
    }
    save(path + ".polylines");
  }

  void compute_barycenter(
    const std::vector<Segment_2>& segments,
    Point_2& b) const {

    FT bx = FT(0), by = FT(0);
    for (const auto& segment : segments) {
      const auto& source = segment.source();
      const auto& target = segment.target();

      bx += source.x(); by += source.y();
      bx += target.x(); by += target.y();
    }
    bx /= static_cast<FT>(segments.size() * 2);
    by /= static_cast<FT>(segments.size() * 2);
    b = Point_2(bx, by);
  }

  void translate_segments(
    const std::vector<Segment_2>& input,
    const Point_2& b,
    std::vector<Segment_2>& segments) const {

    segments.clear();
    segments.reserve(input.size());

    for (const auto& segment : input) {
      const auto& source = segment.source();
      const auto& target = segment.target();

      segments.push_back(Segment_2(
        Point_2(source.x() - b.x(), source.y() - b.y()),
        Point_2(target.x() - b.x(), target.y() - b.y())));
    }
  }

  void compute_bounding_box(
    const std::vector<Segment_2>& segments,
    Point_2& minb, Point_2& maxb) const {

    FT minx = FT(1000000000000), maxx = -FT(1000000000000);
    FT miny = FT(1000000000000), maxy = -FT(1000000000000);
    for (const auto& segment : segments) {
      const auto& source = segment.source();
      const auto& target = segment.target();

      minx = (CGAL::min)(minx, source.x());
      maxx = (CGAL::max)(maxx, source.x());
      miny = (CGAL::min)(miny, source.y());
      maxy = (CGAL::max)(maxy, source.y());

      minx = (CGAL::min)(minx, target.x());
      maxx = (CGAL::max)(maxx, target.x());
      miny = (CGAL::min)(miny, target.y());
      maxy = (CGAL::max)(maxy, target.y());
    }

    const FT d = CGAL::abs(maxx - minx) / FT(10);
    minb = Point_2(minx - d, miny - d);
    maxb = Point_2(maxx + d, maxy + d);
  }

  void set_eps_header(
    const FT llx,
    const FT lly,
    const FT urx,
    const FT ury,
    const std::string title) {

    out << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
    out << "%%BoundingBox: " << llx << " " << lly << " " << urx << " " << ury << std::endl;
    out << "%%Pages: 1" << std::endl;
    out << "%%Creator: Dmitry Anisimov, rudanston@gmail.com" << std::endl;
    out << "%%Title: " << title.c_str() << std::endl;
    out << "%%EndComments" << std::endl;
    out << "%%EndProlog" << std::endl << std::endl;
    out << "%%Page: 1 1" << std::endl << std::endl;
  }

  void add_eps_segment(
    const Segment_2& segment,
    const FT scale,
    const FT line_width,
    const bool dashed) {

    const auto& source = segment.source();
    const auto& target = segment.target();

    out << source.x() * scale << " " << source.y() * scale << " moveto" << std::endl;
    out << target.x() * scale << " " << target.y() * scale << " lineto" << std::endl;
    out << 0 << " setgray" << std::endl;

    if (dashed) out << "[4 1] 0 setdash" << std::endl;
    else out << "[] 0 setdash" << std::endl;

    out << line_width << " setlinewidth" << std::endl;
    out << "stroke" << std::endl << std::endl;
  }

  void add_eps_disc(
    const Point_2& center,
    const FT radius,
    const FT scale) {

    out << 0 << " setgray" << std::endl;
    out << "0 setlinewidth" << std::endl << std::endl;
    out <<
    center.x() * scale << " " <<
    center.y() * scale << " " <<
    radius << " 0 360 arc closepath" << std::endl << std::endl;
    out << "gsave" << std::endl;
    out << 0 << " setgray fill" << std::endl;
    out << "grestore" << std::endl;
    out << "stroke" << std::endl << std::endl;
  }
};

#endif // CGAL_SHAPE_REGULARIZATION_EXAMPLES_SAVER_H
