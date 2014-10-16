// Author(s) : Dmitry Anisimov.
// We use a simple triangle and inexact data type
// in order to test coordinates computed for the center point of the triangle.
// We test the function Compute() and return CGAL::Point_3 type of coordinates.

// It also works with exact kernel.

#include <CGAL/Triangle_coordinates_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef Kernel::Point_3    Point_3;

typedef std::vector<Scalar> Coordinate_vector;

typedef CGAL::Barycentric_coordinates::Tri_coordinates_2<Triangle, Coordinate_vector> Triangle_coordinates;

using std::cout; using std::endl; using std::string;
 
int main()
{
    const Triangle triangle( Point(0, 0), Point(2, 0), Point(1, 2) );

    const Point center(1, 1);

    const Point_3 p = Triangle_coordinates::Compute(triangle, center);

    if( p.x() - (Scalar(1) / Scalar(4)) != Scalar(0) ||
        p.y() - (Scalar(1) / Scalar(4)) != Scalar(0) ||
        p.z() - (Scalar(1) / Scalar(2)) != Scalar(0)  )
    {
        cout << endl << "Triangle_test_with_point_3: FAILED." << endl << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Triangle_test_with_point_3: PASSED." << endl << endl;

    return EXIT_SUCCESS;
}
