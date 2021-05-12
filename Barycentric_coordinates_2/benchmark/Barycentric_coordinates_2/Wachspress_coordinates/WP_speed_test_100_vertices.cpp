// Author: Dmitry Anisimov.
// We test speed of Wachspress coordinates on a set of automatically generated
// points inside a regular polygon with 100 vertices. We use inexact kernel.

#include <CGAL/number_type_config.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;

typedef Scalar_vector::iterator Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Wachspress_2<Kernel> Wachspress;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Wachspress, Kernel> Wachspress_coordinates;

using std::cout; using std::endl; using std::string;

void generate_regular_polygon(const int number_of_vertices, const double polygon_radius, Point_vector &vertices)
{
    const int n = number_of_vertices;
    const double r =  polygon_radius;

    vertices.resize(n);

    for(int i = 0; i < n; ++i)
        vertices[i] = Point(Scalar(r*sin((CGAL_PI / n) + ((i * 2.0 * CGAL_PI) / n))), Scalar(-r*cos((CGAL_PI / n) + ((i * 2.0 * CGAL_PI) / n))));
}

int main()
{
    const int number_of_x_coordinates = 1000;
    const int number_of_y_coordinates = 1000;
    const int number_of_runs          = 1;

    const Scalar one  = Scalar(1);
    const Scalar x_step = one / Scalar(number_of_x_coordinates);
    const Scalar y_step = one / Scalar(number_of_y_coordinates);

    Point_vector vertices;

    const int number_of_vertices = 100;
    const double  polygon_radius = 2;

    generate_regular_polygon(number_of_vertices, polygon_radius, vertices);

    Wachspress_coordinates wachspress_coordinates(vertices.begin(), vertices.end());

    Scalar_vector coordinates(number_of_vertices);
    Overwrite_iterator it = coordinates.begin();

    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) {

        time_to_compute.start();
        for(Scalar x = -one; x <= one; x += x_step) {
            for(Scalar y = -one; y <= one; y += y_step)
                wachspress_coordinates(Point(x, y), it, CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE);
        }
        time_to_compute.stop();

        time += time_to_compute.time();

        time_to_compute.reset();
    }
    const double mean_time = time / number_of_runs;

    cout.precision(10);
    cout << endl << "CPU time to compute Wachspress coordinates (100 vertices) = " << mean_time << " seconds." << endl << endl;

    return EXIT_SUCCESS;
}
