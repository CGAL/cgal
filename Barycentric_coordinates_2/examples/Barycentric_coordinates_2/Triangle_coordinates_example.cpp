// Author(s) : Dmitry Anisimov.

#include <CGAL/Cartesian.h>

#include <CGAL/Triangle_coordinates_2.h>

typedef CGAL::Cartesian<float> Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef std::vector<Scalar>  Coordinate_vector;
typedef std::insert_iterator<Coordinate_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Triangle, Vector_insert_iterator> Triangle_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a triangle.
    const Triangle triangle( Point(0.0f, 0.0f), Point(2.0f, 0.5f), Point(1.0f, 2.0f) );

    // Create an std::vector to store coordinates.
    Coordinate_vector coordinates;

    // Instantiate Triangle coordinates class for the triangle defined above.
    Triangle_coordinates triangle_coordinates(triangle);

    // Print some information about the triangle and coordinate functions.
    triangle_coordinates.print_info();

    // Instantiate some interior, boundary, and exterior points for which we compute coordinates.
    const int number_of_query_points = 18;
    const Point query_points[] = { Point(0.5f , 0.5f ), Point(1.0f, 0.5f ), Point(1.0f , 0.75f), Point(1.0f , 1.0f),                     // interior points
                                   Point(1.0f , 1.25f), Point(1.0f, 1.5f ), Point(0.75f, 1.0f ), Point(1.25f, 1.0f), Point(1.5f, 0.75f),
                                   Point(1.0f , 0.25f), Point(0.5f, 1.0f ), Point(1.5f , 1.25f), Point(1.0f , 2.0f), Point(2.0f, 0.5f ), // boundary points
                                   Point(0.25f, 1.0f ), Point(0.5f, 1.75f), Point(1.5f , 1.75f), Point(1.75f, 1.5f)                      // exterior points
                                 };

    // Reserve memory to store triangle coordinates for 18 query points.
    coordinates.reserve(18 * 3);

    // Compute Triangle coordinates for these points.
    cout << endl << "Computed Triangle coordinates are " << endl << endl;
    for(int i = 0; i < number_of_query_points; ++i) {
        triangle_coordinates.compute(query_points[i], std::inserter(coordinates, coordinates.end()));

        // Output coordinates for each point.
        cout << "Point " << i + 1 << ": ";
        for(int j = 0; j < 3; ++j)
            cout << "coordinate " << j + 1 << " = " << coordinates[i * 3 + j] << "; ";
        cout << endl << endl;
    }

    return EXIT_SUCCESS;
}