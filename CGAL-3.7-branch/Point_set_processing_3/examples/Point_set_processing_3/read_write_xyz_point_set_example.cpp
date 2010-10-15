#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <utility> // defines std::pair
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored as a std::pair.
typedef std::pair<Point, Vector> Pwn;

int main(void)
{
    // Reads a .xyz point set file in points[].
    // Note: read_xyz_points_and_normals() requires an output iterator
    // over points and as well as property maps to access each
    // point position and normal.
    std::vector<Pwn> points;
    std::ifstream in("data/oni.xyz");
    if (!in ||
        !CGAL::read_xyz_points_and_normals(
            in,std::back_inserter(points),
            CGAL::First_of_pair_property_map<Pwn>(),
            CGAL::Second_of_pair_property_map<Pwn>()))
    {
      std::cerr << "Error: cannot read file data/oni.xyz" << std::endl;
      return EXIT_FAILURE;
    }

    // Saves point set.
    // Note: write_xyz_points_and_normals() requires an output iterator
    // over points as well as property maps to access each
    // point position and normal.
    std::ofstream out("oni_copy.xyz");
    if (!out ||
        !CGAL::write_xyz_points_and_normals(
            out, points.begin(), points.end(),
            CGAL::First_of_pair_property_map<Pwn>(),
            CGAL::Second_of_pair_property_map<Pwn>()))
    {
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
