#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Reads a .xyz point set file in points[].
    // Note: read_xyz_point_set() requires an output iterator over points
    //       + property maps to access each point's position and normal.
    //       The position property map has a default value and is omitted here.
    //       The function will skip normal vectors as we do not specify a normal property map.
    std::vector<Point> points;
    std::ifstream in("data/sphere_20k.xyz");
    if (!in ||
        !CGAL::read_xyz_point_set(in,
                                  std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Save point set
    std::ofstream out("sphere_20k_copy.xyz");
    if (!out ||
        !CGAL::write_xyz_point_set(out,
                                   points.begin(), points.end()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

