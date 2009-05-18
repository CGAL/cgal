#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/point_set_property_map.h>
#include <CGAL/IO/read_xyz_points.h>

#include <list>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal; // position + normal vector

int main(void)
{
    // Reads a .xyz point set file in points[].
    std::list<Point_with_normal> points;
    std::ifstream stream("data/sphere_20k.xyz");
    if (!stream ||
        !CGAL::read_xyz_points(stream,
                               std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Estimates normals direction.
    // Note: pca_estimate_normals() requires an iterator over points
    //       + property maps to access each point's position and normal.
    //       The position property map can be omitted here as we use an iterator over Point_3 elements.
    const int nb_neighbors = 7; // K-nearest neighbors
    CGAL::pca_estimate_normals(points.begin(), points.end(),
                               CGAL::make_normal_vector_property_map(points.begin()),
                               nb_neighbors);

    // Orients normals.
    // Note: mst_orient_normals() requires an iterator over points
    //       + property maps to access each point's index, position and normal.
    //       The index and point property maps have default values and are omitted here.
    CGAL::mst_orient_normals(points.begin(), points.end(),
                             CGAL::make_normal_vector_property_map(points.begin()),
                             nb_neighbors);

    return EXIT_SUCCESS;
}

