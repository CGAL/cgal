// Author: Dmitry Anisimov.
// We test speed of triangle coordinates on a set of automatically generated
// points inside a right triangle. We use inexact kernel.

#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;

typedef std::vector<Scalar> Coordinate_vector;
typedef Coordinate_vector::iterator Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> Triangle_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    const int number_of_x_coordinates = 100000;
    const int number_of_y_coordinates = 10000;
    const int number_of_runs          = 1;

    const Scalar zero = Scalar(0);
    const Scalar one  = Scalar(1);
    const Scalar two  = Scalar(2);
    const Scalar x_step = one / Scalar(number_of_x_coordinates);
    const Scalar y_step = one / Scalar(number_of_y_coordinates);

    const Point first_vertex  = Point(zero - x_step, zero - x_step);
    const Point second_vertex = Point(two  + y_step, zero - x_step);
    const Point third_vertex  = Point(zero - x_step, two + y_step );

    Triangle_coordinates triangle_coordinates(first_vertex, second_vertex, third_vertex);

    Coordinate_vector coordinates(3);
    Overwrite_iterator it = coordinates.begin();

    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) {

        time_to_compute.start();
        for(Scalar x = zero; x <= one; x += x_step) {
            for(Scalar y = zero; y <= one; y += y_step)
                triangle_coordinates(Point(x, y), it);
        }
        time_to_compute.stop();

        time += time_to_compute.time();

        time_to_compute.reset();
    }
    const double mean_time = time / number_of_runs;

    cout.precision(10);
    cout << endl << "CPU time to compute Triangle coordinates = " << mean_time << " seconds." << endl << endl;

    return EXIT_SUCCESS;
}
