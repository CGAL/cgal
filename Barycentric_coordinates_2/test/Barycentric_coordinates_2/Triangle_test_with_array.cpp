// Author: Dmitry Anisimov.
// We use a simple triangle and inexact data type
// in order to test coordinates computed for the center point of the triangle.
// We test the function compute_triangle_coordinates_2() and return CGAL::cpp11::array set of coordinates.

// It also works with exact kernel.

#include <cassert>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::cpp11::array<Scalar,3> Triple;

using std::cout; using std::endl; using std::string;
 
int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(2, 0);
    const Point third_vertex  = Point(1, 2);
    const Point center        = Point(1, 1);

    const Triple p = CGAL::Barycentric_coordinates::compute_triangle_coordinates_2(first_vertex, second_vertex, third_vertex, center, Kernel());

    assert(p[0] - (Scalar(1) / Scalar(4)) == Scalar(0) && 
           p[1] - (Scalar(1) / Scalar(4)) == Scalar(0) &&
           p[2] - (Scalar(1) / Scalar(2)) == Scalar(0) );

    if( p[0] - (Scalar(1) / Scalar(4)) != Scalar(0) ||
        p[1] - (Scalar(1) / Scalar(4)) != Scalar(0) ||
        p[2] - (Scalar(1) / Scalar(2)) != Scalar(0)  )
    {
        cout << endl << "Triangle_test_with_array: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Triangle_test_with_array: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
