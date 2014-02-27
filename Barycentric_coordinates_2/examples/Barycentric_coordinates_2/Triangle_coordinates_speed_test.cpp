// Author(s) : Dmitry Anisimov.

#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangle_coordinates_2.h>

// Construct an iterator that takes as an input current data type and pointer to the first element in the standard C++ array.

template<typename Scalar>
    class overwrite_iterator
{
private:
    Scalar* pointer;

public:
    explicit overwrite_iterator(Scalar* _pointer) : pointer(_pointer) { }

    // There are only two operations that we need to overload in order to use the Triangle coordinates class.

    // This operation is intended to return the current coordinate value.
    inline Scalar& operator* () { return *pointer; }

    // This operation is intended to increase the index of the coordinate.
    inline void operator++ () { ++pointer; }
};

// Some convenient typedefs.

typedef CGAL::Real_timer Timer;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT         Scalar;
typedef Kernel::Point_2    Point;
typedef Kernel::Triangle_2 Triangle;

typedef overwrite_iterator<Scalar> Overwrite_iterator;

typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Triangle, Overwrite_iterator> Triangle_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Number of x and y coordinates together gives the number of points.
    const int number_of_x_coordinates = 1000000;
    const int number_of_y_coordinates = 10000;

    // Number of runs to compute the arithmetic mean of the time.
    const int number_of_runs = 10;

    // Compute the uniform step size along x and y directions to change coordinates.
    const Scalar zero = Scalar(0);
    const Scalar one  = Scalar(1);
    const Scalar two  = Scalar(2);
    const Scalar x_step = one / Scalar(number_of_x_coordinates);
    const Scalar y_step = one / Scalar(number_of_y_coordinates);

    // Create a right triangle with a slight offset from zero.
    const Triangle right_triangle( Point(zero - x_step, zero - x_step), Point(two + y_step, zero - x_step), Point(zero - x_step, two + y_step) );

    // Instantiate the Triangle coordinates class for the right triangle defined above.
    Triangle_coordinates triangle_coordinates(right_triangle);

    // Create an instance of a standard C++ array to store coordinates.
    // It has fixed size = 3 = number of vertices.
    Scalar coordinates [3] = {0, 0, 0};

    // Pass pointer to the first element of the array with coordinates in order to overwrite them.
    Overwrite_iterator it( &(coordinates[0]) );

    // Create timer.
    Timer time_to_compute;

    double time = 0.0;
    for(int i = 0; i < number_of_runs; ++i) { // Number of runs

        time_to_compute.start(); // Start clock
        for(Scalar x = zero; x <= one; x += x_step) {
            for(Scalar y = zero; y <= one; y += y_step)
                triangle_coordinates.compute(Point(x, y), it); // Compute 3 coordinate values for each generated point
        }
        time_to_compute.stop(); // Stop clock

        time += time_to_compute.time(); // Add time of the current run to the whole time

        time_to_compute.reset(); // Reset time.
    }

    // Compute the arithmetic mean of all the runs.
    const double mean_time = time / number_of_runs;

    // Output resulting time.
    cout.precision(10);
    cout << endl << "CPU time to compute Triangle coordinates for "
         << number_of_x_coordinates * number_of_y_coordinates << " points = " << mean_time << " seconds.";
    cout << endl << endl;

    return EXIT_SUCCESS;
}