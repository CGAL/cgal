// Author(s) : Dmitry Anisimov.

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Discrete_harmonic_coordinates_2.h>

// Some convenient typedefs.

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::list<Scalar> Coordinate_list;
typedef std::back_insert_iterator<Coordinate_list> List_insert_iterator;
typedef Coordinate_list::const_iterator List_iterator;

typedef CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2<Polygon, List_insert_iterator> Discrete_harmonic_coordinates;
typedef std::pair<List_insert_iterator, bool> Output_type;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a unit square.
    const int number_of_vertices = 4;
    const Point vertices[number_of_vertices] = { Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1) };
    const Polygon unit_square(vertices, vertices + number_of_vertices);

    // Create an std::list to store coordinates.
    Coordinate_list coordinates;

    // Instantiate Discrete Harmonic coordinates class for the unit square defined above.
    Discrete_harmonic_coordinates discrete_harmonic_coordinates(unit_square);

    // Print some information about the polygon and coordinate functions.
    discrete_harmonic_coordinates.print_info();

    // Instantiate the center point of the unit square.
    const Point center(Scalar(1)/Scalar(2), Scalar(1)/Scalar(2));

    // Compute Discrete Harmonic coordinates for the center point.
    // Use parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE.
    Output_type result = discrete_harmonic_coordinates.compute(center, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE);

    // Instantiate other 4 interior points.
    const int number_of_interior_points = 4;
    const Point interior_points[number_of_interior_points] = { Point(Scalar(1)/Scalar(5), Scalar(1)/Scalar(5))  ,
                                                               Point(Scalar(4)/Scalar(5), Scalar(1)/Scalar(5))  ,
                                                               Point(Scalar(4)/Scalar(5), Scalar(4)/Scalar(5))  ,
                                                               Point(Scalar(1)/Scalar(5), Scalar(4)/Scalar(5)) };

    // Compute Discrete Harmonic coordinates for these points and store them at the same list "coordinates" as before.
    for(int i = 0; i < number_of_interior_points; ++i)
        result = discrete_harmonic_coordinates.compute(interior_points[i], std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE);

    // Instantiate 2 boundary points on the second and last edges.
    const Point second_edge(1, Scalar(4)/Scalar(5));
    const Point   last_edge(0, Scalar(4)/Scalar(5));

    // Compute Discrete Harmonic coordinates for these 2 points.
    // Use parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDARY.
    result = discrete_harmonic_coordinates.compute(second_edge, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDARY);
    result = discrete_harmonic_coordinates.compute(last_edge  , std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDARY);

    // Instantiate 2 other boundary points on the first and third edges.
    const Point first_edge(Scalar(1)/Scalar(2), 0);
    const Point third_edge(Scalar(1)/Scalar(2), 1);

    // Compute Discrete Harmonic coordinates using index of an appropriate edge.
    // Do not forget that index counting starts from zero.
    result = discrete_harmonic_coordinates.compute_on_edge(first_edge, 0, std::back_inserter(coordinates));
    result = discrete_harmonic_coordinates.compute_on_edge(third_edge, 2, std::back_inserter(coordinates));

    // Compute Discrete Harmonic coordinates for points at the first and third vertex of the unit square.
    result = discrete_harmonic_coordinates.compute_at_vertex(0, std::back_inserter(coordinates));
    result = discrete_harmonic_coordinates.compute_at_vertex(2, std::back_inserter(coordinates));

    // Instantiate points at the second and fourth vertex of the unit square.
    const Point second_vertex(1, 0);
    const Point fourth_vertex(0, 1);

    // Compute Discrete Harmonic coordinates for these points.
    // Use parameter query_point_location = CGAL::Barycentric_coordinates::AT_VERTEX.
    result = discrete_harmonic_coordinates.compute(second_vertex, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::AT_VERTEX);
    result = discrete_harmonic_coordinates.compute(fourth_vertex, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::AT_VERTEX);

    // Instantiate 2 points outside the unit square - one from the left and one from the right.
    const Point left_most(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(2));
    const Point right_most(Scalar(3)/Scalar(2), Scalar(1)/Scalar(2));

    // Compute Discrete Harmonic coordinates for these 2 points.
    // Use parameter query_point_location = CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE.
    result = discrete_harmonic_coordinates.compute(left_most , std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE);
    result = discrete_harmonic_coordinates.compute(right_most, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE);

    // Output computed coordinates.
    cout << endl << "Exact Discrete Harmonic coordinates for all the defined points: " << endl << endl;
    int index; List_iterator it;
    for(index = 0, it = coordinates.begin(); it != coordinates.end(); ++index, ++it) {
        cout << "Coordinate " << index % number_of_vertices + 1 << " = " << *it << " ";
        if((index +  1) %      number_of_vertices  == 0) cout << endl;
        if((index + 13) % (4 * number_of_vertices) == 0) cout << endl;
    }

    // Return status of the last computation.
    const string status = (result.second == true ? "SUCCESS." : "FAILURE.");
    cout << endl << "Status of the last computation: " << status << endl << endl;

    return EXIT_SUCCESS;
}