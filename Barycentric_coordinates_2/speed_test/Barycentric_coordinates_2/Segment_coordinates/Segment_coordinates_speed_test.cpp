// Author(s) : Dmitry Anisimov.
// We test speed of Segment coordinates on a set of automatically generated
// points inside a unit segment. We use inexact kernel.

#include <CGAL/Real_timer.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Segment_coordinates_2.h>

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT        Scalar;
typedef Kernel::Point_2   Point;
typedef Kernel::Segment_2 Segment;

typedef std::vector<Scalar> Coordinate_vector;
typedef Coordinate_vector::iterator Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Segment_coordinates_2<Segment, Overwrite_iterator> Segment_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    const int number_of_iterations = 1000000;
    const int number_of_points     = 10000;
    const int number_of_runs       = 1;

    const Scalar zero = Scalar(0);
    const Scalar one  = Scalar(1);
    const Scalar step = one / Scalar(number_of_points);

    const Segment unit_segment( Point(zero - step, zero), Point(one + step, zero) );

    Segment_coordinates segment_coordinates(unit_segment);

    Coordinate_vector coordinates;
    coordinates.resize(2);
    Overwrite_iterator it = coordinates.begin();

    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) {

        time_to_compute.start();
        for(int j = 0; j < number_of_iterations; ++j) {
            for(Scalar x = zero; x <= one; x += step) 
                segment_coordinates.compute(Point(x, zero), it);
        }
        time_to_compute.stop();

        time += time_to_compute.time();

        time_to_compute.reset();
    }
    const double mean_time = time / number_of_runs;

    cout.precision(10);
    cout << endl << "CPU time to compute Segment coordinates = " << mean_time << " seconds." << endl << endl;
    
    return EXIT_SUCCESS;
}
