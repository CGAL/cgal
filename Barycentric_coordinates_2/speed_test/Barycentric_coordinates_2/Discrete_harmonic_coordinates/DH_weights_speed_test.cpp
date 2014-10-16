// Author(s) : Dmitry Anisimov.
// We test speed of Discrete Harmonic weights on a set of automatically generated
// points inside a unit square. We use inexact kernel.

#include <CGAL/Real_timer.h>

#include <CGAL/Polygon_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Discrete_harmonic_coordinates_2.h>

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::Polygon_2<Kernel> Polygon;

typedef std::vector<Scalar> Weight_vector;
typedef Weight_vector::iterator Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2<Polygon, Overwrite_iterator> Discrete_harmonic_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    const int number_of_x_coordinates = 1000000;
    const int number_of_y_coordinates = 1000;
    const int number_of_runs          = 1;

    const Scalar zero = Scalar(0);
    const Scalar one  = Scalar(1);
    const Scalar x_step = one / Scalar(number_of_x_coordinates);
    const Scalar y_step = one / Scalar(number_of_y_coordinates);

    const Point vertices[4] = { Point(zero - x_step, zero - y_step),
                                Point(one  + x_step, zero - y_step),
                                Point(one  + x_step, one  + y_step),
                                Point(zero - x_step, one  + y_step)
                              };
    const Polygon unit_square(vertices, vertices + 4);

    Discrete_harmonic_coordinates discrete_harmonic_coordinates(unit_square);

    Weight_vector weights;
    weights.resize(4);
    Overwrite_iterator it = weights.begin();

    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) {

        time_to_compute.start();
        for(Scalar x = zero; x <= one; x += x_step) {
            for(Scalar y = zero; y <= one; y += y_step)
                discrete_harmonic_coordinates.compute_weights(Point(x, y), it);
        }
        time_to_compute.stop();

        time += time_to_compute.time();

        time_to_compute.reset();
    }
    const double mean_time = time / number_of_runs;

    cout.precision(10);
    cout << endl << "CPU time to compute Discrete Harmonic weights (4 vertices) = " << mean_time << " seconds." << endl << endl;

    return EXIT_SUCCESS;
}
