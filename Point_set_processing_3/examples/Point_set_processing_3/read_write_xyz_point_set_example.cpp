#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Reads a .xyz point set file in points[].
    // Note: read_xyz_points() requires an output iterator over points
    //       + a property map to access each point's position.
    //       The position property map can be omitted here as we use an iterator over Point_3 elements.
    std::vector<Point> points;
    std::ifstream in("data/sphere_20k.xyz");
    if (!in ||
        !CGAL::read_xyz_points(in,
                               std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Saves point set.
    // Note: write_xyz_points() requires an output iterator over points
    //       + a property map to access each point's position.
    //       The position property map can be omitted here as we use an iterator over Point_3 elements.
    std::ofstream out("sphere_20k_copy.xyz");
    if (!out ||
        !CGAL::write_xyz_points(out,
                                points.begin(), points.end()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

