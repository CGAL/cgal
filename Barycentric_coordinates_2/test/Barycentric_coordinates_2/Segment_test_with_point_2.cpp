// Author(s) : Dmitry Anisimov.
// We use a simple segment of length 2 and inexact data type
// in order to test coordinates computed for the center point of the segment.
// We test the function Compute() and return CGAL::Point_2 type of coordinates.

// It also works with exact kernel.

#include <CGAL/Segment_coordinates_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT        Scalar;
typedef Kernel::Point_2   Point;
typedef Kernel::Segment_2 Segment;

typedef std::vector<Scalar> Coordinate_vector;

typedef CGAL::Barycentric_coordinates::Seg_coordinates_2<Segment, Coordinate_vector> Segment_coordinates;

using std::cout; using std::endl; using std::string;
 
int main()
{
    const Segment segment( Point(0, 0), Point(2, 0) );

    const Point center(1, 0);

    const Point p = Segment_coordinates::Compute(segment, center);

    const Scalar half = Scalar(1)/Scalar(2);
    if( p.x() - half != Scalar(0) ||
        p.y() - half != Scalar(0)  )
    {
        cout << endl << "Segment_test_with_point_2: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Segment_test_with_point_2: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
