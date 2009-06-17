#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair; 

int main(void)
{
    // Reads a .xyz point set file in points[].
    // Note: read_xyz_points_and_normals() requires an output iterator over points
    // + property maps to access each point's position and normal.
    std::vector<PointVectorPair> points;
    std::ifstream in("data/oni.xyz");
    if (!in ||
        !CGAL::read_xyz_points_and_normals(
            in,std::back_inserter(points),
            CGAL::First_of_pair_property_map<PointVectorPair>(),
            CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
      return EXIT_FAILURE;
    }

    // Saves point set.
    // Note: write_xyz_points_and_normals() requires an output iterator over points
    // + property maps to access each point's position and normal.
    std::ofstream out("oni_copy.xyz");
    if (!out ||
        !CGAL::write_xyz_points_and_normals(
            out, points.begin(), points.end(),
            CGAL::First_of_pair_property_map<PointVectorPair>(),
            CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

