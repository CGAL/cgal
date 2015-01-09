#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

// Some convenient typedefs.
typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Scalar_vector;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a triangle.
    const Point  first_vertex(0.0f, 0.0f);
    const Point second_vertex(2.0f, 0.5f);
    const Point  third_vertex(1.0f, 2.0f);

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;

    // Instantiate the class Triangle_coordinates_2 for the triangle defined above.
    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    // Print some information about the triangle and coordinates.
    triangle_coordinates.print_information();

    // Instantiate some interior, boundary, and exterior query points for which we compute coordinates.
    const int number_of_query_points = 18;
    const Point query_points[] = { Point(0.5f , 0.5f ), Point(1.0f, 0.5f ), Point(1.0f , 0.75f), Point(1.0f , 1.0f),                     // interior query points
                                   Point(1.0f , 1.25f), Point(1.0f, 1.5f ), Point(0.75f, 1.0f ), Point(1.25f, 1.0f), Point(1.5f, 0.75f),
                                   Point(1.0f , 0.25f), Point(0.5f, 1.0f ), Point(1.5f , 1.25f), Point(1.0f , 2.0f), Point(2.0f, 0.5f ), // boundary query points
                                   Point(0.25f, 1.0f ), Point(0.5f, 1.75f), Point(1.5f , 1.75f), Point(1.75f, 1.5f)                      // exterior query points
                                 };

    // Reserve memory to store triangle coordinates for 18 query points.
    coordinates.reserve(number_of_query_points * 3);

    // Compute triangle coordinates for these points.
    cout << endl << "Computed triangle coordinates: " << endl << endl;
    for(int i = 0; i < number_of_query_points; ++i) {

        triangle_coordinates(query_points[i], std::inserter(coordinates, coordinates.end()));

        // Output the coordinates for each point.
        cout << "Point " << i + 1 << ": ";
        for(int j = 0; j < 3; ++j)
            cout << "coordinate " << j + 1 << " = " << coordinates[i * 3 + j] << "; ";
        cout << endl << endl;
    }

    return EXIT_SUCCESS;
}
