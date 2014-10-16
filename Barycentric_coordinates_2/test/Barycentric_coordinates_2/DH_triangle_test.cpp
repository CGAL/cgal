// Author(s) : Dmitry Anisimov.
// In this test we compute Discrete Harmonic coordinates for ~2400 strictly interior points
// with respect to a triangle and compare them with those from Triangle coordinates.
// They must be the same. The chosen data type is exact.

// Does not work with inexact kernel. We get inconsistency when comparing Triangle and Discrete Harmonic coordinates.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Triangle_coordinates_2.h>
#include <CGAL/Discrete_harmonic_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Tri_coordinates_2<Triangle, Coordinate_vector>        Triangle_coordinates;
typedef CGAL::Barycentric_coordinates::DH_coordinates_2<Polygon, Coordinate_vector> Discrete_harmonic_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Triangle tri_triangle( Point(0, 0), Point(1, 0), Point(0, 1) );

    Triangle_coordinates triangle_coordinates(tri_triangle);

    const Point vertices[3] = { Point(0, 0), Point(1, 0), Point(0, 1) };

    const Polygon dh_triangle(vertices, vertices + 3);

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(dh_triangle);

    Coordinate_vector tri_coordinates;
    Coordinate_vector  dh_coordinates;

    const Scalar step  = Scalar(1) / Scalar(100);
    const Scalar scale = Scalar(50);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type tri_result = triangle_coordinates.compute(point, tri_coordinates);
            const Output_type  dh_result = discrete_harmonic_coordinates.compute(point, dh_coordinates);

            if( tri_coordinates[count + 0] - dh_coordinates[count + 0] != Scalar(0) ||
                tri_coordinates[count + 1] - dh_coordinates[count + 1] != Scalar(0) ||
                tri_coordinates[count + 2] - dh_coordinates[count + 2] != Scalar(0)  )
            {
                cout << endl << tri_result.second << " " << dh_result.second << " DH_triangle_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "DH_triangle_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
