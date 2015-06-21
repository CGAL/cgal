// Author: Dmitry Anisimov.
// In this test we compute discrete harmonic coordinates for ~2400 strictly interior points
// with respect to a triangle and compare them with those from triangle coordinates.
// They must be the same. The chosen data type is exact.

// Does not work with inexact kernel. We get inconsistency when comparing triangle and discrete harmonic coordinates.

#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::vector<Point>  Point_vector;

typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;
typedef CGAL::Barycentric_coordinates::Discrete_harmonic_2<Kernel> Discrete_harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Discrete_harmonic, Kernel> Discrete_harmonic_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Point first_vertex  = Point(0, 0);
    const Point second_vertex = Point(1, 0);
    const Point third_vertex  = Point(0, 1);

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    Point_vector vertices(3);
    vertices[0] = first_vertex; vertices[1] = second_vertex; vertices[2] = third_vertex;

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertices.begin(), vertices.end());

    Coordinate_vector tri_coordinates;
    Coordinate_vector  dh_coordinates;

    const Scalar step  = Scalar(1) / Scalar(100);
    const Scalar scale = Scalar(50);

    int count = 0;
    const Scalar limit = scale*step;

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type tri_result = triangle_coordinates(point, tri_coordinates);
            const Output_type  dh_result = discrete_harmonic_coordinates(point, dh_coordinates);

            assert(tri_coordinates[count + 0] - dh_coordinates[count + 0] == Scalar(0) &&
                   tri_coordinates[count + 1] - dh_coordinates[count + 1] == Scalar(0) && 
                   tri_coordinates[count + 2] - dh_coordinates[count + 2] == Scalar(0) );

            if( tri_coordinates[count + 0] - dh_coordinates[count + 0] != Scalar(0) ||
                tri_coordinates[count + 1] - dh_coordinates[count + 1] != Scalar(0) ||
                tri_coordinates[count + 2] - dh_coordinates[count + 2] != Scalar(0)  )
            {
                cout << endl << "DH_triangle_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "DH_triangle_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
