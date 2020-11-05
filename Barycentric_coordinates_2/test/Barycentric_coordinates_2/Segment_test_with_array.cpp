// Author: Dmitry Anisimov.
// We use a simple segment of length 2 and inexact data type
// in order to test coordinates computed for the center point of the segment.
// We test the function compute_segment_coordinates_2() and return std::array set of coordinates.

// It also works with exact kernel.

#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::array<Scalar,2> Pair;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(2, 0);
    const Point center        = Point(1, 0);

    const Pair p = CGAL::Barycentric_coordinates::compute_segment_coordinates_2(first_vertex, second_vertex, center, Kernel());

    const Scalar half = Scalar(1)/Scalar(2);

    assert( p[0] - half == Scalar(0) && p[1] - half == Scalar(0) );

    if( p[0] - half != Scalar(0) ||
        p[1] - half != Scalar(0)  )
    {
        cout << endl << "Segment_test_with_array: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Segment_test_with_array: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
