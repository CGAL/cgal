#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>

#include <list>
#include <iostream>
#include <fstream>
#include <cstring>

using std::strcmp;

// leda_rational, or Gmpq, or Quotient<MP_float>
using NT = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<NT>;
using Segment_traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Traits_2 = CGAL::Arr_curve_data_traits_2<Segment_traits_2, int>;
using Point_2 = Traits_2::Point_2;
using Segment_2 = Segment_traits_2::Curve_2;
using Curve_2 = Traits_2::Curve_2;
using Diagram_1 = CGAL::Envelope_diagram_1<Traits_2>;
using Curve_list = std::list<Curve_2>;

enum Coord_input_format {
  F_RATIONAL,            // Coordinates given as rational numbers.
  F_INTEGER,             // Coordinates given as integers.
  F_DOUBLE               // Coordinates given as doubles (with decimal points).
};

/*! reads a set of line segments from an input file.
 * \param filename The name of the input file.
 * \param format The coordinate format.
 * \param segs Output: The segments.
 * \return Whether the segments were successfully read.
 */
bool read_segments(const char* filename, Coord_input_format format, Curve_list& segs) {
  segs.clear();

  // Open the input file.
  std::ifstream ifile(filename);

  if (! ifile.is_open()) {
    std::cerr << "Failed to open <" << filename << ">." << std::endl;
    return (false);
  }

  // Read the segments.
  unsigned int n_segments;
  ifile >> n_segments;

  unsigned int k;
  for (k = 0; k < n_segments; k++) {
    NT x1, y1, x2, y2;

    // Read the coordinates of the current segment.
    if (format == F_INTEGER) {
      int ix1, iy1, ix2, iy2;
      ifile >> ix1 >> iy1 >> ix2 >> iy2;
      x1 = NT(ix1);
      y1 = NT(iy1);
      x2 = NT(ix2);
      y2 = NT(iy2);
    }
    else if (format == F_DOUBLE) {
      const int denom = 100000;
      double dx1, dy1, dx2, dy2;
      ifile >> dx1 >> dy1 >> dx2 >> dy2;
      x1 = NT(static_cast<int>(denom * dx1), denom);
      y1 = NT(static_cast<int>(denom * dy1), denom);
      x2 = NT(static_cast<int>(denom * dx2), denom);
      y2 = NT(static_cast<int>(denom * dy2), denom);
    }
    else {
      ifile >> x1 >> y1 >> x2 >> y2;
    }

    Segment_2 seg = Segment_2(Point_2(x1, y1), Point_2(x2, y2));
    segs.push_back(Curve_2(seg, static_cast<int>(segs.size())));
  }
  ifile.close();

  return (true);
}

/*! checks if a curve with the same associated data as the input
 * curve is in the given range.
 * \param begin The beginning of the range.
 * \param end Past-the-end iterator.
 * \param c The curve, the data of which we are searching.
 * \return True if we found a curve with the same data.
 */
template <typename I>
bool find_curve(I begin, I end, const Curve_2& c) {
  while (begin != end) {
    if (begin->data() == c.data()) return true;
    ++begin;
  }
  return false;
}

/*! checks the envelope of a given set of segments.
 * \param segs The segments.
 * \param diag The diagram.
 * \param is_lower Does the diagram represent the lower or the upper envelope.
 * \return Whether the diagram structure is correct.
 */
bool check_envelope(const Curve_list& segs, const Diagram_1& diag, bool is_lower) {
  // Go over minimization diagram.
  Diagram_1::Edge_const_handle e = diag.leftmost();
  Diagram_1::Vertex_const_handle v;

  // If the diagram is empty, the segment set must also be empty.
  if (e == diag.rightmost()) return (segs.empty());

  // Start from the first finite edge.
  Kernel ker;
  Kernel::Construct_midpoint_2 midpoint = ker.construct_midpoint_2_object();
  Kernel::Construct_min_vertex_2 min_ver = ker.construct_min_vertex_2_object();
  Kernel::Construct_max_vertex_2 max_ver = ker.construct_max_vertex_2_object();
  Kernel::Compare_y_at_x_2 comp_y_at_x = ker.compare_y_at_x_2_object();
  Point_2 p_mid;
  CGAL::Comparison_result res1, res2;
  CGAL::Comparison_result y_res;
  Curve_list::const_iterator sit;

  v = diag.right_vertex(e);
  e = diag.right_edge(v);

  while (e != diag.rightmost()) {
    const Point_2& p_left = diag.point(diag.left_vertex(e));
    const Point_2& p_right = diag.point(diag.right_vertex(e));

    // Get the midpoint of the current edge.
    p_mid = midpoint(p_left, p_right);

    const auto& e_curves = diag.edge_curves(e);
    const auto& v_curves = diag.vertex_curves(v);

    // Check that all associated curves are equal at p_mid:
    for (auto it = e_curves.begin(); it != e_curves.end(); ++it) {
      if (comp_y_at_x(p_mid, e_curves.front(), *it) != CGAL::EQUAL) {
        std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with "
                  << "both of the following curves which are not equal: "
                  << " [" << e_curves.front() << "], [" << *it << "]" << std::endl;
        return false;
      }
    }

    // Check that all associated curves are equal at the vertex:
    const Point_2& p_v = diag.point(v);
    for (auto it = v_curves.begin(); it != v_curves.end(); ++it) {
      if (comp_y_at_x(p_v, *it) != CGAL::EQUAL) {
        std::cerr << "The vertex (" << p_v << ") is associated with "
                  << " [" << *it << "] but they are not equal at the vertex." << std::endl;
        return false;
      }
    }

    // Go over all segments.
    for (sit = segs.begin(); sit != segs.end(); ++sit) {
      // Check if p_mid lies in the x-range of the current segment.
      res1 = CGAL::compare_x(p_mid, min_ver(*sit));
      res2 = CGAL::compare_x(p_mid, max_ver(*sit));

      if (res1 != CGAL::SMALLER && res2 != CGAL::LARGER) {
        // Check the diagram edge.
        if (diag.empty_edge_curves(e)) {
          // If p_mid is in the x-range of the given segment, the current edge
          // cannot be empty.
          std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is empty, "
                    << "but the segment [" << *sit << "] is defined over it." << std::endl;
          return (false);
        }
        else {
          // Compare the segment with the segment associated with the edge.
          y_res = comp_y_at_x(p_mid, *sit, e_curves.front());

          // Check for violations of the diagram properties.
          if ((is_lower && y_res == CGAL::SMALLER) || (! is_lower && y_res == CGAL::LARGER)) {
            std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with the segment ["
                      << e_curves.front() << "], but the segment [" << *sit << "] violates the envelope properties."
                      << std::endl;
            return (false);
          }

          // If it is equal, the segment should be in the diagram.
          if (y_res == CGAL::EQUAL && find_curve(e_curves.begin(), e_curves.end(), *sit) == false) {
            std::cerr << "The edge (" << p_left << ") -> (" << p_right << ") is associated with the segment ["
                      << e_curves.front() << "], but the segment [" << *sit << "] is not in the list of its curves."
                      << std::endl;
            return (false);
          }
        }
      }

      // Check if v lies in the x-range of the current segment.
      res1 = CGAL::compare_x(p_v, min_ver(*sit));
      res2 = CGAL::compare_x(p_v, max_ver(*sit));

      if (res1 != CGAL::SMALLER && res2 != CGAL::LARGER) {
        // Check the diagram vertex.
        if (v_curves.empty()) {
          // If v is in the x-range of the given segment, the current vertex
          // cannot be empty.
          std::cerr << "The vertex (" << p_v << ") is empty, " << "but the segment [" << *sit << "] is defined over it."
                    << std::endl;
          return (false);
        }
        else {
          // Compare the segment with the segment associated with the vertex.
          y_res = comp_y_at_x(p_v, *sit);

          // Check for violations of the diagram properties.
          if ((is_lower && y_res == CGAL::LARGER) || (! is_lower && y_res == CGAL::SMALLER)) {
            std::cerr << "The vertex (" << p_v << ") and the segment [" << *sit << "] violate the envelope properties."
                      << std::endl;
            return (false);
          }

          // If it is equal, the segment should be in the diagram.
          if (y_res == CGAL::EQUAL && find_curve(v_curves.begin(), v_curves.end(), *sit) == false) {
            std::cerr << "The vertex (" << p_v << ") does not contain the segment ["
                      << *sit << "] in the list of its curves, but they have" << " equal values." << std::endl;
            return (false);
          }
        }
      }

    }
    v = diag.right_vertex(e);
    e = diag.right_edge(v);
  }

  // If we reached here, the diagram is valid.
  return (true);
}

/*!
 * The main program.
 */
int main(int argc, char *argv[]) {
  if (argc % 2 == 0 || argc == 1) {
    std::cerr << "Usage: " << argv[0] << " <input file> [ -q | -i | -d ]" << std::endl;
    return (1);
  }

  int number_of_tests = (argc - 1) / 2;
  for (int i = 0; i < number_of_tests; ++i) {
    std::cout << "Checking input file: " << argv[2*i + 1] << std::endl;

    // Determine the input format.
    Coord_input_format   format = F_RATIONAL;

    if (strcmp(argv[2*i + 2], "-i") == 0 || strcmp(argv[2*i + 2], "-I") == 0) format = F_INTEGER;
    else if (strcmp(argv[2*i + 2], "-d") == 0 || strcmp(argv[2*i + 2], "-D") == 0) format = F_DOUBLE;

    // Read the input segments.
    Curve_list segments;
    if (! read_segments(argv[2*i + 1], format, segments)) return (1);

    // Compute their lower envelope.
    Diagram_1 min_diag;

    CGAL::lower_envelope_2(segments.begin(), segments.end(), min_diag);

    // Check the lower envelope.
    if (! check_envelope(segments, min_diag, true)) {
      std::cerr << "Problems in the lower-envelope computation." << std::endl;
      return (1);
    }
    else {
      std::cout << "The lower envelope is valid." << std::endl;
    }

    // Compute the upper envelope.
    Diagram_1              max_diag;

    CGAL::upper_envelope_2(segments.begin(), segments.end(), max_diag);

    // Check the upper envelope.
    if (! check_envelope (segments, max_diag, false)) {
      std::cerr << "Problems in the upper-envelope computation." << std::endl;
      return (1);
    }
    else {
      std::cout << "The upper envelope is valid." << std::endl;
    }
  }
  return (0);
}
