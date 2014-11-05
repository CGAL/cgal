#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <utility> // defines std::pair
#include <vector>
#include <cstdlib>
#include <fstream>

#include <CGAL/Voronoi_covariance_3/vcm_estimate_edges.h>

#define RR 0.1
#define rr 0
#define TT 0.2

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

int main (int argc, char *argv[]) {
    PointList points;
    if (argc < 5) {
        std::cerr << "Not enough arguments" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1];

    // Read the off file passed in argument
    std::cout << "Loading " << input_filename << std::endl;
    std::ifstream stream(input_filename.c_str());
    bool success = false;
    success = stream &&
        CGAL::read_off_points(stream,
                              std::back_inserter(points),
                              CGAL::First_of_pair_property_map<PointVectorPair>());
    if (!success) {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }

    // Sharp edges estimation via VCM
    double R, r, threshold;
    R = atof(argv[2]);
    r = atof(argv[3]);
    threshold = atof(argv[4]);
    std::cout << "Edges estimation via VCM" << std::endl;
    std::vector<Segment> polylines;
    polylines = vcm_estimate_edges(points.begin(), points.end(),
                                   CGAL::First_of_pair_property_map<PointVectorPair>(),
                                   R, r, threshold,
                                   Kernel());

    std::ofstream file_polylines("polylines");
    for (std::vector<Segment>::iterator it = polylines.begin();
         it != polylines.end();
         ++it) {
        file_polylines << *it << "\n";
    }

    return 0;
}

