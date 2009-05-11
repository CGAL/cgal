#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/random_simplify_point_set.h>
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

    // Randomly simplifies.
    const double removed_percentage = 75.0; // removed percentage
    CGAL::random_simplify_point_set(points.begin(), points.end(),
                                    removed_percentage);

    return EXIT_SUCCESS;
}

