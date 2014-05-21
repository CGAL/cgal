#include <CGAL/Barycentric_traits_2.h>
#include <CGAL/Barycentric_coordinates_2.h>
#include <CGAL/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// Namespace alias.
namespace BC = CGAL::Barycentric_coordinates;

// Some convenient typedefs.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef BC::Barycentric_traits_2<Kernel> Barycentric_traits;

typedef Barycentric_traits::FT      Scalar;
typedef Barycentric_traits::Point_2 Point;

typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;

typedef Point_vector::iterator InputIterator;
typedef std::back_insert_iterator<Scalar_vector> Vector_insert_iterator;
typedef std::pair<Vector_insert_iterator, bool> Output_type;

typedef BC::Discrete_harmonic_coordinates_2<Barycentric_traits> Discrete_harmonic;
typedef BC::Barycentric_coordinates_2<InputIterator, Discrete_harmonic, Barycentric_traits> Discrete_harmonic_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a unit square.
    const int number_of_vertices = 4;
    Point_vector vertex;
    vertex.resize(number_of_vertices);
    vertex[0] = Point(0, 0); vertex[1] = Point(1, 0); vertex[2] = Point(1, 1); vertex[3] = Point(0, 1);

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;

    // Instantiate the class Discrete_harmonic_coordinates_2 for the unit square defined above.
    Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertex.begin(), vertex.end());

    // Print some information about the polygon and coordinate functions.
    discrete_harmonic_coordinates.print_information();

    // Instantiate the center point of the unit square.
    const Point center(Scalar(1)/Scalar(2), Scalar(1)/Scalar(2));

    // Compute discrete harmonic coordinates for the center point.
    // Use the parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE.
    Output_type result = discrete_harmonic_coordinates.compute(center, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE);

    // Instantiate other 4 interior points.
    const int number_of_interior_points = 4;
    const Point interior_points[number_of_interior_points] = { Point(Scalar(1)/Scalar(5), Scalar(1)/Scalar(5))  ,
                                                               Point(Scalar(4)/Scalar(5), Scalar(1)/Scalar(5))  ,
                                                               Point(Scalar(4)/Scalar(5), Scalar(4)/Scalar(5))  ,
                                                               Point(Scalar(1)/Scalar(5), Scalar(4)/Scalar(5)) };

    // Compute discrete harmonic coordinates for these points and store them at the same vector "coordinates" as before.
    for(int i = 0; i < number_of_interior_points; ++i)
        result = discrete_harmonic_coordinates.compute(interior_points[i], std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE);

    // Instantiate 2 boundary points on the second and last edges.
    const Point second_edge(1, Scalar(4)/Scalar(5));
    const Point   last_edge(0, Scalar(4)/Scalar(5));

    // Compute discrete harmonic coordinates for these 2 points.
    // Use the parameter query_point_location = CGAL::Barycentric_coordinates::ON_BOUNDARY.
    result = discrete_harmonic_coordinates.compute(second_edge, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDARY);
    result = discrete_harmonic_coordinates.compute(last_edge  , std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_BOUNDARY);

    // Instantiate 2 other boundary points on the first and third edges.
    const Point first_edge(Scalar(1)/Scalar(2), 0);
    const Point third_edge(Scalar(1)/Scalar(2), 1);

    // Compute discrete harmonic coordinates using index of an appropriate edge.
    // Do not forget that index counting starts from zero.
    result = discrete_harmonic_coordinates.compute_on_edge(first_edge, 0, std::back_inserter(coordinates));
    result = discrete_harmonic_coordinates.compute_on_edge(third_edge, 2, std::back_inserter(coordinates));

    // Compute discrete harmonic coordinates for the points at the first and third vertex of the unit square.
    result = discrete_harmonic_coordinates.compute_on_vertex(0, std::back_inserter(coordinates));
    result = discrete_harmonic_coordinates.compute_on_vertex(2, std::back_inserter(coordinates));

    // Instantiate points at the second and fourth vertex of the unit square.
    const Point second_vertex(1, 0);
    const Point fourth_vertex(0, 1);

    // Compute discrete harmonic coordinates for these points.
    // Use the parameter query_point_location = CGAL::Barycentric_coordinates::ON_VERTEX.
    result = discrete_harmonic_coordinates.compute(second_vertex, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_VERTEX);
    result = discrete_harmonic_coordinates.compute(fourth_vertex, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_VERTEX);

    // Instantiate 2 points outside the unit square - one from the left and one from the right.
    const Point left_most(Scalar(-1)/Scalar(2), Scalar(1)/Scalar(2));
    const Point right_most(Scalar(3)/Scalar(2), Scalar(1)/Scalar(2));

    // Compute discrete harmonic coordinates for these 2 points.
    // Use the parameter query_point_location = CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE.
    result = discrete_harmonic_coordinates.compute(left_most , std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE);
    result = discrete_harmonic_coordinates.compute(right_most, std::back_inserter(coordinates), CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE);

    // Output the computed coordinate values.
    cout << endl << "Exact discrete harmonic coordinates for all the defined points: " << endl << endl;
    const size_t number_of_query_points = coordinates.size(); 
    for(int index = 0; index < int(number_of_query_points); ++index) {
        cout << "Coordinate " << index % number_of_vertices + 1 << " = " << coordinates[index] << " ";
        if((index +  1) %      number_of_vertices  == 0) cout << endl;
        if((index + 13) % (4 * number_of_vertices) == 0) cout << endl;
    }

    // Return status of the last computation.
    const string status = (result.second == true ? "SUCCESS." : "FAILURE.");
    cout << endl << "Status of the last computation: " << status << endl << endl;

    return EXIT_SUCCESS;
}