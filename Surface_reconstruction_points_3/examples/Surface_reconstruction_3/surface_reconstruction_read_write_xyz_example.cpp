#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>
#include <CGAL/IO/surface_reconstruction_write_xyz.h>

#include <deque>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::deque<Point_with_normal> PointList;

int main(void)
{
    // Read a .xyz point set file in points[]
    PointList points;
    if (!CGAL::surface_reconstruction_read_xyz("data/sphere_20k.xyz",
                                               std::back_inserter(points),
                                               false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Save point set
    if(!CGAL::surface_reconstruction_write_xyz("sphere_20k_copy.xyz",
                                               points.begin(), points.end()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

