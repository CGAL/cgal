#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/pca_normal_estimation.h>
#include <CGAL/mst_normal_orientation.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <deque>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Orientable_normal_3<Kernel> Orientable_normal; // normal vector + orientation
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal; // position + normal vector
typedef std::deque<Point_with_normal> PointList;

int main(void)
{
    // Read a .xyz point set file in points[]
    PointList points;
    std::ifstream stream("data/sphere_20k.xyz");
    if (!stream || 
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points),
                                  false /*skip normals*/))
    {
      return EXIT_FAILURE;
    }

    // Estimate normals direction.
    std::deque<Orientable_normal> output;
    const int nb_neighbors = 7; // K-nearest neighbors
    CGAL::pca_normal_estimation(points.begin(), points.end(),
                                std::back_inserter(output),
                                nb_neighbors);

    // Orient normals.
    // mst_normal_orientation() requires an iterator over points
    // + property maps to access each point's index, position and normal.
    // We use the points index as iterator.
    boost::identity_property_map index_id; // identity
    CGAL::mst_normal_orientation(
           (std::size_t)0, points.size(), // use the points index as iterator
           index_id, // index -> index property map = identity
           boost::make_iterator_property_map(points.begin(), index_id), // index -> position prop. map
           boost::make_iterator_property_map(output.begin(), index_id), // index -> normal prop. map
           nb_neighbors);

    return EXIT_SUCCESS;
}

