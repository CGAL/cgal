// Author(s) : Dmitry Anisimov.
// In this test we compute Mean Value coordinates for ~2400 strictly interior points with respect to a triangle
// and compare them with those from Triangle coordinates. They must be the same. The chosen
// data type is inexact and epsilon = 1.0e-14.

// Works with an exact type, too.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangle_coordinates_2.h>
#include <CGAL/Mean_value_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Coordinate_vector;
typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Tri_coordinates_2<Triangle, Coordinate_vector> Triangle_coordinates;
typedef CGAL::Barycentric_coordinates::MV_coordinates_2<Polygon, Coordinate_vector> Mean_value_coordinates;

typedef std::pair<Vector_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    const Triangle tri_triangle( Point(0, 0), Point(1, 0), Point(0, 1) );

    Triangle_coordinates triangle_coordinates(tri_triangle);

    const Point vertices[3] = { Point(0, 0), Point(1, 0), Point(0, 1) };

    const Polygon mv_triangle(vertices, vertices + 3);

    Mean_value_coordinates mean_value_coordinates(mv_triangle);

    Coordinate_vector tri_coordinates;
    Coordinate_vector  mv_coordinates;

    const Scalar step  = Scalar(1) / Scalar(100);
    const Scalar scale = Scalar(50);

    int count = 0;
    const Scalar limit = scale*step;
    const Scalar epsilon = Scalar(1) / Scalar(std::pow(10.0, 14.0));

    for(Scalar x = step; x < limit; x += step) {
        for(Scalar y = step; y < limit; y += step) {
            const Point point(x, y);

            const Output_type tri_result = triangle_coordinates.compute(point, tri_coordinates);
            const Output_type  mv_result = mean_value_coordinates.compute(point, mv_coordinates);

            if( (tri_coordinates[count + 0] - mv_coordinates[count + 0]) > epsilon ||
                (tri_coordinates[count + 1] - mv_coordinates[count + 1]) > epsilon ||
                (tri_coordinates[count + 2] - mv_coordinates[count + 2]) > epsilon  )
            {
                cout << endl << tri_result.second << " " << mv_result.second << " MV_triangle_test: FAILED." << endl << endl;
                exit(EXIT_FAILURE);
            }
            count += 3;
        }
    }

    cout << endl << "MV_triangle_test: PASSED." << endl << endl;
    
    return EXIT_SUCCESS;
}
