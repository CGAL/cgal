#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

using Kernel    = CGAL::Simple_cartesian<double>;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;

int main() {

  // Create input segments.
  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0.2, 0.0), Point_2(1.2, 0.0)),
    Segment_2(Point_2(1.2, 0.1), Point_2(2.2, 0.1)),
    Segment_2(Point_2(2.2, 0.0), Point_2(2.0, 2.0)),
    Segment_2(Point_2(2.0, 2.0), Point_2(1.0, 2.0)),
    Segment_2(Point_2(1.0, 1.9), Point_2(0.0, 1.9)),
    Segment_2(Point_2(0.0, 2.0), Point_2(0.2, 0.0))
  };

  // Regularize all segments: both angles and offsets.
  CGAL::Shape_regularization::Segments::
    regularize_segments(segments);
}
