#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <deque>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

int main(void)
{
    // Reads a .xyz point set file 
    std::deque<Point> points;
    std::ifstream stream("data/sphere_20k.xyz");
    if (!stream || 
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points),
                                  false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Computes average spacing
    const unsigned int nb_neighbors = 7;
    typedef std::deque<Point>::iterator Iterator;
    FT average_spacing = CGAL::compute_average_spacing<Iterator>(points.begin(), points.end(),
                                                                 nb_neighbors);
    std::cout << "Average spacing: " << average_spacing << std::endl;

    return EXIT_SUCCESS;
}

