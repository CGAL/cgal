#include <CGAL/Barycentric_traits_2.h>
#include <CGAL/Segment_coordinates_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// Namespace alias.
namespace BC = CGAL::Barycentric_coordinates;

// Some convenient typedefs.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef BC::Barycentric_traits_2<Kernel> Barycentric_traits;

typedef Barycentric_traits::FT      Scalar;
typedef Barycentric_traits::Point_2 Point;

typedef BC::Segment_coordinates_2<Barycentric_traits> Segment_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a segment.
    const Point vertex_1 = Point(0, Scalar(2)/Scalar(5));
    const Point vertex_2 = Point(2, Scalar(2)/Scalar(5));

    // Instantiate three interior and two exterior query points.
    const Point query_points[5] = { Point(Scalar(2) /Scalar(5), Scalar(2)/Scalar(5)), // interior query points
                                    Point(1                   , Scalar(2)/Scalar(5)),
                                    Point(Scalar(8) /Scalar(5), Scalar(2)/Scalar(5)),
                                    Point(Scalar(-1)/Scalar(5), Scalar(2)/Scalar(5)), // exterior query points
                                    Point(Scalar(11)/Scalar(5), Scalar(2)/Scalar(5))
                                  };

    // Compute segment coordinates for all the defined points.
    // We use a static function and return segment coordinates stored in the CGAL::Point_2 data structure.
    cout << endl << "Computed segment coordinates are " << endl << endl;
    for(int i = 0; i < 5; ++i) {
        const Point point = Segment_coordinates::static_compute(vertex_1, vertex_2, query_points[i], Barycentric_traits());

        // Output both coordinates for each point.
        cout << "Point " << i + 1 << " = (" << point.x() << ", " << point.y() << ");" << endl;
    }
    cout << endl;

    return EXIT_SUCCESS;
}