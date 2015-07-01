#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>

#include <utility> // defines std::pair
#include <vector>
#include <fstream>

#include <boost/foreach.hpp>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

typedef CGAL::cpp11::array<double,6> Covariance;

int main (int , char**) {
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

    // Estimates covariance matrices per points.
    double R = 0.2,
           r = 0.1;
    std::vector<Covariance> cov;
    CGAL::First_of_pair_property_map<PointVectorPair> point_pmap;

    CGAL::compute_vcm(points.begin(), points.end(), point_pmap, cov, R, r, Kernel());

    // Find the points on the edges.
    // Note that this step is not expensive and can be done several time to get better results
    double threshold = 0.16;
    std::ofstream output("points_on_edges.xyz");
    int i = 0;
    BOOST_FOREACH(const PointVectorPair& p, points)
    {
      if (CGAL::vcm_is_on_feature_edge(cov[i], threshold))
          output << p.first << "\n";
      ++i;
    }

    return 0;
}

