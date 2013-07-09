#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/denoise_point_set.h>
#include <CGAL/Timer.h>
//#include <CGAL/Memory_sizer.h>

#include <utility> // defines std::pair
#include <list>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

int main(void)
{
    // Reads a .xyz point set file in points[].
    std::list<PointVectorPair> points;
    std::ifstream stream("data/sphere_20k_normal.xyz");
    //std::ifstream stream("data/oni.xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(stream,
                       std::back_inserter(points),
                       CGAL::First_of_pair_property_map<PointVectorPair>(),
                       CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
        std::cerr << "Error: cannot read file data/sphere_20k_normal.xyz" << std::endl;
        return EXIT_FAILURE;
    }

    // 
    int k = 200;
    int iter_number = 5;

    CGAL::Timer task_timer;
    task_timer.start();
    std::cout << "Run algorithm example: " << std::endl;

    for (int i = 0; i < iter_number; i++)
    {
      CGAL::denoise_points_with_normals(points.begin(), points.end(),
            CGAL::First_of_pair_property_map<PointVectorPair>(),
            CGAL::Second_of_pair_property_map<PointVectorPair>(),
            k);
    }

    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, " ;
   /*   << (memory>>20) << " Mb allocated" << std::endl;*/
    task_timer.stop();  



    //// Saves point set.
    //// Note: write_xyz_points_and_normals() requires an output iterator
    //// over points as well as property maps to access each
    //// point position and normal.
    std::ofstream out("data/sphere_20k_denoised.xyz");  
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

