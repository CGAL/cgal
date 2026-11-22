#include <fstream>
#include <cstring>

#include "test_utils.h"

using std::strcmp;

enum Coord_input_format {
  F_RATIONAL,            // coordinates given as rational numbers.
  F_INTEGER,             // coordinates given as integers.
  F_DOUBLE               // coordinates given as doubles (with decimal points).
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
  std::size_t n_segments;
  ifile >> n_segments;

  for (std::size_t k = 0; k < n_segments; ++k) {
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

/*!
 * The main program.
 */
int main(int argc, char* argv[]) {
  if (argc % 2 == 0 || argc == 1) {
    std::cerr << "Usage: " << argv[0] << " <input file> [ -q | -i | -d ]" << std::endl;
    return (1);
  }

  int number_of_tests = (argc - 1) / 2;
  for (int i = 0; i < number_of_tests; ++i) {
    std::cout << "Checking input file: " << argv[2*i + 1] << std::endl;

    // Determine the input format.
    Coord_input_format format = F_RATIONAL;

    if (strcmp(argv[2*i + 2], "-i") == 0 || strcmp(argv[2*i + 2], "-I") == 0) format = F_INTEGER;
    else if (strcmp(argv[2*i + 2], "-d") == 0 || strcmp(argv[2*i + 2], "-D") == 0) format = F_DOUBLE;

    // Read the input segments.
    Curve_list segments;
    if (! read_segments(argv[2*i + 1], format, segments)) return (1);

    // Compute their lower envelope.
    Diagram_1 min_diag;

    CGAL::Envelope_2::lower_envelope_2(segments.begin(), segments.end(), min_diag);

    // Check the lower envelope.
    if (! check_envelope(segments, min_diag, true)) {
      std::cerr << "Problems in the lower-envelope computation." << std::endl;
      return (1);
    }
    else {
      std::cout << "The lower envelope is valid." << std::endl;
    }

    // Compute the upper envelope.
    Diagram_1 max_diag;

    CGAL::Envelope_2::upper_envelope_2(segments.begin(), segments.end(), max_diag);

    // Check the upper envelope.
    if (! check_envelope(segments, max_diag, false)) {
      std::cerr << "Problems in the upper-envelope computation." << std::endl;
      return (1);
    }
    else {
      std::cout << "The upper envelope is valid." << std::endl;
    }
  }
  return (0);
}
