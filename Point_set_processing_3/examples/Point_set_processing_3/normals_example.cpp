#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>

#include <utility> // defines std::pair
#include <list>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/sphere_1k.xyz";
    // Reads a .xyz point set file in points[].
    std::list<PointVectorPair> points;
    std::ifstream stream(fname);
    if (!stream ||
        !CGAL::read_xyz_points(stream,
                               std::back_inserter(points),
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())))
    {
      std::cerr << "Error: cannot read file " << fname<< std::endl;
        return EXIT_FAILURE;
    }

    // Estimates normals direction.
    // Note: pca_estimate_normals() requiresa range of points
    // as well as property maps to access each point's position and normal.
    const int nb_neighbors = 18; // K-nearest neighbors = 3 rings
    CGAL::pca_estimate_normals<Concurrency_tag>
      (points, nb_neighbors,
       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
       normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    // Orients normals.
    // Note: mst_orient_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
    std::list<PointVectorPair>::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(points, nb_neighbors,
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                               normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    // Optional: delete points with an unoriented normal
    // if you plan to call a reconstruction algorithm that expects oriented normals.
    points.erase(unoriented_points_begin, points.end());

    return EXIT_SUCCESS;
}

