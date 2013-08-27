#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/tags.h>

#include <utility> // defines std::pair
#include <list>
#include <string>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

int main(void)
{
    const std::string INPUT_FILENAME_WITHOUT_EXT = "data/qtr_piston_noise_with_normal";

    // Reads a .xyz point set file in points[].
    std::vector<PointVectorPair> points;
    std::ifstream stream(INPUT_FILENAME_WITHOUT_EXT + ".xyz");
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(stream,
                       std::back_inserter(points),
                       CGAL::First_of_pair_property_map<PointVectorPair>(),
                       CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
       std::cerr << "Error: cannot read file " 
                 << INPUT_FILENAME_WITHOUT_EXT << ".xyz" << std::endl;
       return EXIT_FAILURE;
    }

    //Algorithm parameters
    int k = 120;  //neighborhood size
    double sharpness_sigma = 25; // control sharpness(0-90), 
                                 // the bigger, the result will be smoother.
    int iter_number = 3; //times of projection

    CGAL::Timer task_timer;
    task_timer.start();
    std::cout << "Run algorithm example: " << std::endl;

    for (int i = 0; i < iter_number; i++)
    {
      std::cout << std::endl << "Iteration: " << i+1 << std::endl;
     
      double error = 
      CGAL::bilateral_smooth_point_set <CGAL::Parallel_tag>(
            points.begin(), 
            points.end(),
            CGAL::First_of_pair_property_map<PointVectorPair>(),
            CGAL::Second_of_pair_property_map<PointVectorPair>(),
            k,
            sharpness_sigma);

      std::cout << std::endl << "move error: " << error << std::endl;
    }

    long memory = CGAL::Memory_sizer().virtual_size();
    std::cout << "done: " << task_timer.time() << " seconds, "
              << (memory>>20) << " Mb allocated" << std::endl;
    task_timer.stop();  


    //// Save point set.
    std::ofstream out(INPUT_FILENAME_WITHOUT_EXT + "_jet DENOISED.xyz");  
    std::ofstream out("data/sphere_20k_denoised.xyz");  
    if (!out ||
      !CGAL::write_xyz_points_and_normals(
      out, points.begin(), points.end(), 
      CGAL::First_of_pair_property_map<PointVectorPair>(),
      CGAL::Second_of_pair_property_map<PointVectorPair>()))
    {
      return EXIT_FAILURE;
    }

    system("Pause");

    return EXIT_SUCCESS;
}

