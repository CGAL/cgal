#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/average_spacing_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

#include <deque>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Read a .xyz point set file in points[]
    std::deque<Point> points;
    if (!CGAL::surface_reconstruction_read_xyz("data/sphere_20k.xyz",
                                               std::back_inserter(points),
                                               false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Average Spacing
    const unsigned int nb_neighbors = 7;
    typedef std::deque<Point>::iterator Iterator;
    FT average_spacing = CGAL::average_spacing_3<Iterator,FT>(points.begin(), points.end(),
                                                              nb_neighbors);

    return EXIT_SUCCESS;
}

