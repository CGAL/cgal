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
    // Reads a .xyz point set file
    std::deque<Point> points;
    std::ifstream stream("data/oni.xyz");
    if (!stream || 
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points),
                                  false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Randomly simplifies
    const double removed_percentage = 75.0; // removed percentage
    CGAL::random_simplify_point_set(
                    points.begin(), points.end(),
                    removed_percentage);

    return EXIT_SUCCESS;
}

