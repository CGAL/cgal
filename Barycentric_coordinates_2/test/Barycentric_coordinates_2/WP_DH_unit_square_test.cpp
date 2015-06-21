// Author: Dmitry Anisimov.
// In this test we compute Wachspress and discrete harmonic coordinates for ~9800 strictly interior
// points from the unit square and compare values of both coordinate functions.
// Since the vertices of this polygon lie on a circle, the coordinate values must be the same.
// The chosen data type is exact.

// Does not work with inexact kernel. We get inconsistency when comparing coordinates.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Wachspress_2<Kernel> Wachspress;
typedef CGAL::Barycentric_coordinates::Discrete_harmonic_2<Kernel> Discrete_harmonic;

typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Wachspress, Kernel> Wachspress_coordinates;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Discrete_harmonic, Kernel> Discrete_harmonic_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    Point_vector vertices(4);

    vertices[0] = Point(0, 0); vertices[1] = Point(1, 0);
    vertices[2] = Point(1, 1); vertices[3] = Point(0, 1);

    Wachspress_coordinates               wachspress_coordinates(vertices.begin(), vertices.end());
    Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertices.begin(), vertices.end());

    Coordinate_vector wp_coordinates;
    Coordinate_vector dh_coordinates;

    const Scalar step  = Scalar(1) / Scalar(100);
    const Scalar scale = Scalar(100);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type wp_result = wachspress_coordinates(point, wp_coordinates);
            const Output_type dh_result = discrete_harmonic_coordinates(point, dh_coordinates);

            assert(wp_coordinates[count + 0] - dh_coordinates[count + 0] == Scalar(0) &&
                   wp_coordinates[count + 1] - dh_coordinates[count + 1] == Scalar(0) &&
                   wp_coordinates[count + 2] - dh_coordinates[count + 2] == Scalar(0) );

            if( wp_coordinates[count + 0] - dh_coordinates[count + 0] != Scalar(0) ||
                wp_coordinates[count + 1] - dh_coordinates[count + 1] != Scalar(0) ||
                wp_coordinates[count + 2] - dh_coordinates[count + 2] != Scalar(0)  )
            {
                cout << endl << "WP_DH_unit_square_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "WP_DH_unit_square_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
