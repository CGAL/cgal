#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>

#include <utility> // defines std::pair
#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

int main (int argc, char *argv[]) {
    // Reads a .xyz point set file in points[].
    std::list<PointVectorPair> points;
    std::ifstream stream("data/fandisk.off");
    if (!stream ||
        !CGAL::read_off_points(stream,
                               std::back_inserter(points),
                               CGAL::First_of_pair_property_map<PointVectorPair>()))
    {
        std::cerr << "Error: cannot read file data/fandisk.off" << std::endl;
        return EXIT_FAILURE;
    }

    // Estimates feature edges.
    double R = 0.2,
           r = 0.1,
           threshold = 0.16;
    std::vector<Segment> polylines;
    polylines = vcm_estimate_edges(points.begin(), points.end(),
                                   CGAL::First_of_pair_property_map<PointVectorPair>(),
                                   R, r, threshold,
                                   Kernel());

    return 0;
}

