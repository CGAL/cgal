#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/outlier_removal_3.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <deque>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Read a .xyz point set file in points[]
    std::deque<Point> points;
    std::ifstream stream("data/oni.xyz");
    if (!stream || 
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points),
                                  false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Outlier Removal
    const double threshold_percent_avg_knn_sq_dst = 5.0 /* % */; // percentage of outliers to remove
    const int nb_neighbors = 7; // K-nearest neighbors
    std::deque<Point> output;
    CGAL::outlier_removal_3(
                    points.begin(), points.end(),
                    std::back_inserter(output),
                    nb_neighbors,
                    threshold_percent_avg_knn_sq_dst);

    return EXIT_SUCCESS;
}

