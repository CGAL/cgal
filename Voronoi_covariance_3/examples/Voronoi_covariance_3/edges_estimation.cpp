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

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

// VCM
typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

int main (int argc, char *argv[]) {
    PointList points;
    if (argc < 3) {
        std::cerr << "Not enough arguments" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

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
    R = 0.1;
    r = 0;
    threshold = 0.5;
    std::cout << "Edges estimation via VCM" << std::endl;
    vcm_estimate_edges(points.begin(), points.end(),
                       CGAL::First_of_pair_property_map<PointVectorPair>(),
                       R, r, threshold,
                       Kernel(),
                       Covariance());


    // Save the point set with normals
    std::cout << "Writing " << output_filename << std::endl;
    std::ofstream os(output_filename.c_str());
    if (!os ||
        !CGAL::write_xyz_points_and_normals(os,
                                            points.begin(), points.end(),
                                            CGAL::First_of_pair_property_map<PointVectorPair>(),
                                            CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return 1;
    }

    return 0;
}

