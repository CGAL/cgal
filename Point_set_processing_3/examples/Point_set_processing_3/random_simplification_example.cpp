#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/random_simplify_point_set.h>
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

    // Random Point Set Simplification options
    const double random_simplification_percentage = 50.0 /* % */; // percentage of outliers to remove
    std::deque<Point> output;
    CGAL::random_simplify_point_set(
                    points.begin(), points.end(),
                    std::back_inserter(output),
                    random_simplification_percentage);

    return EXIT_SUCCESS;
}

