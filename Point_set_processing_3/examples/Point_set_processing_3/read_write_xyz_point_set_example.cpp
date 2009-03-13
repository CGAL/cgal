#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include <deque>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::deque<Point_with_normal> PointList;

int main(void)
{
    // Read a .xyz point set file in points[]
    PointList points;
    if (!CGAL::read_xyz_point_set("data/sphere_20k.xyz",
                                  std::back_inserter(points),
                                  false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Save point set
    if(!CGAL::write_xyz_point_set("sphere_20k_copy.xyz",
                                  points.begin(), points.end()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

