// Author(s) : Dmitry Anisimov.
// In this test we compute Wachspress and Discrete Harmonic coordinates for ~9800 strictly interior
// points from the unit square and compare values of both coordinate functions.
// Since the vertices of this polygon lie on a circle, the coordinate values must be the same.
// The chosen data type is exact.

// Does not work with inexact kernel. We get inconsistency when comparing coordinates.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Wachspress_coordinates_2.h>
#include <CGAL/Discrete_harmonic_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::WP_coordinates_2<Polygon, Coordinate_vector> Wachspress_coordinates;
typedef CGAL::Barycentric_coordinates::DH_coordinates_2<Polygon, Coordinate_vector> Discrete_harmonic_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point vertices[4] = { Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1) };

    const Polygon unit_square(vertices, vertices + 4);

    Wachspress_coordinates               wachspress_coordinates(unit_square);
    Discrete_harmonic_coordinates discrete_harmonic_coordinates(unit_square);

    Coordinate_vector wp_coordinates;
    Coordinate_vector dh_coordinates;

    const Scalar step  = Scalar(1) / Scalar(100);
    const Scalar scale = Scalar(100);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type wp_result = wachspress_coordinates.compute(point, wp_coordinates);
            const Output_type dh_result = discrete_harmonic_coordinates.compute(point, dh_coordinates);

            if( wp_coordinates[count + 0] - dh_coordinates[count + 0] != Scalar(0) ||
                wp_coordinates[count + 1] - dh_coordinates[count + 1] != Scalar(0) ||
                wp_coordinates[count + 2] - dh_coordinates[count + 2] != Scalar(0)  )
            {
                cout << endl << wp_result.second << " " << dh_result.second << " WP_DH_unit_square_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "WP_DH_unit_square_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
