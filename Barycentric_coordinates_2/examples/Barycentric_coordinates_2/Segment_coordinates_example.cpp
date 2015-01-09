#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

// Namespace alias.
namespace BC = CGAL::Barycentric_coordinates;

// Some convenient typedefs.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;

typedef CGAL::cpp11::array<Scalar,2> Pair;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a segment.
    const Point  first_vertex(0, Scalar(2)/Scalar(5));
    const Point second_vertex(2, Scalar(2)/Scalar(5));

    // Instantiate three interior and two exterior query points.
    const Point query_points[5] = { Point(Scalar(2) /Scalar(5), Scalar(2)/Scalar(5)), // interior query points
                                    Point(1                   , Scalar(2)/Scalar(5)),
                                    Point(Scalar(8) /Scalar(5), Scalar(2)/Scalar(5)),
                                    Point(Scalar(-1)/Scalar(5), Scalar(2)/Scalar(5)), // exterior query points
                                    Point(Scalar(11)/Scalar(5), Scalar(2)/Scalar(5))
                                  };

    // Compute segment coordinates for all the defined points.
    // We use a global function and return the segment coordinates stored in an array of the type CGAL::cpp11::array<FT,2>.
    cout << endl << "Computed segment coordinates: " << endl << endl;
    for(int i = 0; i < 5; ++i) {
        const Pair pair = BC::compute_segment_coordinates_2(first_vertex, second_vertex, query_points[i], Kernel());

        // Output both coordinates for each point.
        cout << "Pair of coordinates # " << i + 1 << " = (" << pair[0] << ", " << pair[1] << ");" << endl;
    }
    cout << endl;

    return EXIT_SUCCESS;
}
