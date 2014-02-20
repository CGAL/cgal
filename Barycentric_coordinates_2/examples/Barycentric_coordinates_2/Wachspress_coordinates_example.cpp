// Author(s) : Dmitry Anisimov.

#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Random.h>
#include <CGAL/convex_hull_2.h>

#include <CGAL/Wachspress_coordinates_2.h>

// Some convenient typedefs.

typedef CGAL::Cartesian<double> Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;

typedef CGAL::Barycentric_coordinates::WP_coordinates_2<Polygon> Wachspress_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Choose how many random points we want to generate.
    const int number_of_points = 1000;

    // Choose a diameter of generated set of points.
    const double diameter = 1.0;

    // Create a random numbers generator.
    CGAL::Random rand;

    // Create vectors to store x and y coordinates of generated points.
    Scalar_vector x, y;

    x.resize(number_of_points);
    y.resize(number_of_points);

    // Create vectors to store generated points and vertices of a convex polygon.
    Point_vector points, vertices;

    points.resize(number_of_points);

    // Generate a set of random points between 0 and diameter.
    for(int i = 0; i < number_of_points; ++i) {

        x[i] = rand.get_double(0.0, diameter);
        y[i] = rand.get_double(0.0, diameter);

        points[i] = Point(Scalar(x[i]), Scalar(y[i]));
    }

    // Find a convex hull of the generated set of points.
    // This convex hull gives us vertices of a convex polygon containing all the generated points.
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(vertices));

    const int number_of_vertices = vertices.size();

    // Construct a convex polygon.
    const Polygon convex_polygon(vertices.begin(), vertices.end());

    // Instantiate Wachspress coordinates class for the convex polygon defined above.
    Wachspress_coordinates wachspress_coordinates(convex_polygon);

    // Print some information about the polygon and coordinate functions.
    wachspress_coordinates.print_info();
    
    // Compute Wachspress coordinates for all the defined points.
    cout << endl << "Computed Wachspress coordinates are " << endl << endl;
    for(int i = 0; i < number_of_points; ++i) {
        // Compute coordinates.
        Scalar_vector coordinates;
        coordinates.reserve(number_of_vertices);
        wachspress_coordinates.compute(points[i], coordinates);

        // Output computed coordinates.
        cout << "Point " << i + 1 << ": " << endl;
        for(int j = 0; j < number_of_vertices; ++j)  cout << "Coordinate " << j + 1 << " = " << coordinates[j] << "; " << endl;
        cout << endl;
    }

    return EXIT_SUCCESS;
}