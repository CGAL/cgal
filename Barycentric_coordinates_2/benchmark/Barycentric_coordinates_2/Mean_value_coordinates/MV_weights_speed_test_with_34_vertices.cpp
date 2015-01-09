// Author: Dmitry Anisimov.
// We test speed of mean value weights on a set of automatically generated
// points inside a concave polygon with 34 vertices. We use inexact kernel.

#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef std::vector<Scalar> Weight_vector;
typedef std::vector<Point>  Point_vector;

typedef Weight_vector::iterator Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Mean_value_2<Kernel> Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel> Mean_value_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    const int number_of_x_coordinates = 100000;
    const int number_of_y_coordinates = 1000;
    const int number_of_runs          = 1;

    const Scalar zero = Scalar(0);
    const Scalar one  = Scalar(1);
    const Scalar x_step = one / Scalar(number_of_x_coordinates);
    const Scalar y_step = one / Scalar(number_of_y_coordinates);

    Point_vector vertices(34);

    vertices[0]  = Point(zero, zero - y_step);                              vertices[1] = Point(one, zero - y_step);
    vertices[2]  = Point(Scalar(3)    / Scalar(2), Scalar(-1) / Scalar(2)); vertices[3] = Point(2, Scalar(-1) / Scalar(2));
    vertices[4]  = Point(Scalar(5)    / Scalar(2), 0);                      vertices[5] = Point(2, Scalar(1)  / Scalar(2));
    vertices[6]  = Point(Scalar(5)    / Scalar(2), 1);                      vertices[7] = Point(3, Scalar(3)  / Scalar(4));
    vertices[8]  = Point(3, Scalar(5) / Scalar(4));                         vertices[9] = Point(Scalar(5)   / Scalar(2), Scalar(7)  / Scalar(4));
    vertices[10] = Point(3, Scalar(5) / Scalar(2));                         vertices[11] = Point(Scalar(5)  / Scalar(2), Scalar(5)  / Scalar(2));
    vertices[12] = Point(Scalar(9)    / Scalar(4), 2);                      vertices[13] = Point(Scalar(7)  / Scalar(4), 2);
    vertices[14] = Point(2, Scalar(5) / Scalar(2));                         vertices[15] = Point(Scalar(3)  / Scalar(2), Scalar(5)  / Scalar(2));
    vertices[16] = Point(Scalar(5)    / Scalar(4), 2);                      vertices[17] = Point(Scalar(3)  / Scalar(4), 2);
    vertices[18] = Point(1, Scalar(5) / Scalar(2));                         vertices[19] = Point(Scalar(1)  / Scalar(2), Scalar(5)  / Scalar(2));
    vertices[20] = Point(Scalar(1)    / Scalar(4), 2);                      vertices[21] = Point(Scalar(-1) / Scalar(4), 2);
    vertices[22] = Point(zero, Scalar(5)  / Scalar(2));                     vertices[23] = Point(Scalar(-1) / Scalar(2), Scalar(5)  / Scalar(2));
    vertices[24] = Point(Scalar(-3)   / Scalar(4), 2);                      vertices[25] = Point(Scalar(-1) / Scalar(2), Scalar(3)  / Scalar(2));
    vertices[26] = Point(Scalar(-5)   / Scalar(4), Scalar(3)  / Scalar(2)); vertices[27] = Point(Scalar(-1) / Scalar(2), Scalar(3)  / Scalar(4));
    vertices[28] = Point(-one, Scalar(1)  / Scalar(2));                     vertices[29] = Point(-one, zero);
    vertices[30] = Point(Scalar(-3)   / Scalar(2), zero);                   vertices[31] = Point(Scalar(-3) / Scalar(2), Scalar(-1) / Scalar(2));
    vertices[32] = Point(Scalar(-1)   / Scalar(2), Scalar(-1) / Scalar(2)); vertices[33] = Point(Scalar(-1) / Scalar(2), zero - y_step);

    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    Weight_vector weights(34);
    Overwrite_iterator it = weights.begin();

    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) {

        time_to_compute.start();
        for(Scalar x = zero; x <= one; x += x_step) {
            for(Scalar y = zero; y <= one; y += y_step)
                mean_value_coordinates.compute_weights(Point(x, y), it);
        }
        time_to_compute.stop();

        time += time_to_compute.time();

        time_to_compute.reset();
    }
    const double mean_time = time / number_of_runs;

    cout.precision(10);
    cout << endl << "CPU time to compute Mean Value weights (34 vertices) = " << mean_time << " seconds." << endl << endl;
    
    return EXIT_SUCCESS;
}
