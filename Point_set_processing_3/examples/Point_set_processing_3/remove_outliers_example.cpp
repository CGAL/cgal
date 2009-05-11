#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Reads a .xyz point set file in points[].
    std::vector<Point> points;
    std::ifstream stream("data/oni.xyz");
    if (!stream ||
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Removes outliers
    const double removed_percentage = 5.0; // removed percentage
    const int nb_neighbors = 7; // considers 7 nearest neighbor points
    CGAL::remove_outliers(points.begin(), points.end(),
                          nb_neighbors,
                          removed_percentage);

    return EXIT_SUCCESS;
}

