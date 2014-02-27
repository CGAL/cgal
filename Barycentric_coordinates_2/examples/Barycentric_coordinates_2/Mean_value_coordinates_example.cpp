// Author(s) : Dmitry Anisimov.

#include <deque>

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mean_value_coordinates_2.h>

// Some convenient typedefs.

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::deque<Scalar> Deque_container;
typedef std::front_insert_iterator<Deque_container> Deque_insert_iterator;

typedef CGAL::Barycentric_coordinates::Mean_value_coordinates_2<Polygon, Deque_insert_iterator> Mean_value_coordinates;
typedef std::pair<Deque_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{    
    // Construct a star-shaped polygon.
    const int number_of_vertices = 10;
    const Point vertices[] = { Point(0.0, 0.0), Point(0.1, -0.8), Point(0.3, 0.0), Point(0.6, -0.5), Point(0.6 , 0.1)  ,
                               Point(1.1, 0.6), Point(0.3,  0.2), Point(0.1, 0.8), Point(0.1,  0.2), Point(-0.7, 0.0) };
    const Polygon star_shaped_polygon(vertices, vertices + number_of_vertices);

    // Create an std::deque to store coordinates.
    Deque_container coordinates;

    // Instantiate Mean Value coordinates class for the polygon defined above.
    Mean_value_coordinates mean_value_coordinates(star_shaped_polygon);

    // Print some information about the polygon and coordinate functions.
    mean_value_coordinates.print_info();

    // Instantiate some interior points in the polygon.
    const int number_of_interior_points = 8;
    const Point interior_points[] = { Point(0.12, -0.45), Point(0.55, -0.3), Point(0.9 , 0.45),
                                      Point(0.15,  0.35), Point(-0.4, 0.04), Point(0.11, 0.11),
                                      Point(0.28,  0.12), // the only point in the kernel of the star_shaped_polygon
                                      Point(0.55,  0.11)
                                    };

    // Compute Mean Value coordinates for all the defined interior points.

    // We speed up the computation using the O(n) algorithm called with parameter CGAL::Barycentric_coordinates::FAST.
    const CGAL::Barycentric_coordinates::Type_of_algorithm type_of_algorithm = CGAL::Barycentric_coordinates::FAST;

    // Use parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE.
    const CGAL::Barycentric_coordinates::Query_point_location query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE;

    for(int i = 0; i < number_of_interior_points; ++i) {
        const Output_type result = mean_value_coordinates.compute(interior_points[i], std::front_inserter(coordinates), query_point_location, type_of_algorithm);

        // Output coordinates for each point.
        const string status = (result.second == true ? "SUCCESS." : "FAILURE.");
        cout << endl << "For the point " << i + 1 << " status of the computation: " << status << endl;

        for(int j = 0; j < number_of_vertices; ++j)
            cout << "Coordinate " << j + 1 << " = " << coordinates[number_of_vertices - j - 1] << endl;
    }

    // If we need only the unnormalized weights for the last point, we can compute them as follows.

    // Instantiate an std::deque to store weights.
    Deque_container weights;

    // Compute weights.
    const int last_point_index = number_of_interior_points - 1;
    const Output_type result = mean_value_coordinates.compute_weights(interior_points[last_point_index], std::front_inserter(weights));

    // Compute their sum.
    Scalar mv_denominator = Scalar(0);
    for(int j = 0; j < number_of_vertices; ++j) mv_denominator += weights[j];

    // Invert this sum.
    const Scalar mv_inverted_denominator = Scalar(1) / mv_denominator;

    // Output weights.
    const string status = (result.second == true ? "SUCCESS." : "FAILURE.");
    cout << endl << "Status of the weights' computation for the point " << last_point_index + 1 << ": " << status << endl;

    for(int j = 0; j < number_of_vertices; ++j)
        cout << "Weight " << j + 1 << " = " << weights[number_of_vertices - j - 1] << endl;

    // Now, if we normalize the weight functions, we recover values of the Mean Value coordinates for the last point computed earlier.
    cout << endl << "After normalization, for the point " << last_point_index + 1 << " Mean Value coordinates are " << endl;
    for(int j = 0; j < number_of_vertices; ++j)
        cout << "Coordinate " << j + 1 << " = " << weights[number_of_vertices - j - 1] * mv_inverted_denominator << endl;
    cout << endl;

    return EXIT_SUCCESS;
}