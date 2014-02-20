// Author(s) : Dmitry Anisimov.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Segment_coordinates_2.h>

// Some convenient typedefs.

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT        Scalar;
typedef Kernel::Point_2   Point;
typedef Kernel::Segment_2 Segment;

typedef CGAL::Barycentric_coordinates::Seg_coordinates_2<Segment> Segment_coordinates;

using std::cout; using std::endl; using std::string;

int main()
{
    // Construct a segment.
    const Segment segment( Point(0, Scalar(2)/Scalar(5)), Point(2, Scalar(2)/Scalar(5)) );

    // Instantiate three interior and two exterior points.
    const Point query_points[5] = { Point(Scalar(2) /Scalar(5), Scalar(2)/Scalar(5)), // interior points
                                    Point(1                   , Scalar(2)/Scalar(5)),
                                    Point(Scalar(8) /Scalar(5), Scalar(2)/Scalar(5)),
                                    Point(Scalar(-1)/Scalar(5), Scalar(2)/Scalar(5)), // exterior points
                                    Point(Scalar(11)/Scalar(5), Scalar(2)/Scalar(5))
                                  };

    // Compute Segment coordinates for all the defined points.
    // We use a static function and return Segment coordinates stored in the CGAL::Point_2 data structure.
    cout << endl << "Computed Segment coordinates are " << endl << endl;
    for(int i = 0; i < 5; ++i) {
        const Point point = Segment_coordinates::Compute(segment, query_points[i]);

        // Output both coordinates for each point.
        cout << "For point  " << i + 1 << ": "   << endl << endl;
        cout << "Coordinate " << 1 << " = " << point.x() << "; " << endl;
        cout << "Coordinate " << 2 << " = " << point.y() << "; " << endl;
        cout << endl;
    }

    return EXIT_SUCCESS;
}