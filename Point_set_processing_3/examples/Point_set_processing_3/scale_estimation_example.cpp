#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;

int main (int argc, char** argv)
{

  const char* fname = (argc>1)?argv[1]:"data/sphere_20k.xyz";
  
  CGAL::Timer task_timer; 
  
  std::ifstream stream(fname);
  if (stream)
    {
      std::vector<Point_3> points;
      std::cerr << "TRYING 3D CASE" << std::endl;
      if (CGAL::read_xyz_points(stream, std::back_inserter(points)))
        {
          task_timer.start();
          
          // estimate global scale
          std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale (points.begin(), points.end());
          FT range_scale = CGAL::estimate_global_range_scale (points.begin(), points.end());

          std::size_t memory = CGAL::Memory_sizer().virtual_size();
          double time = task_timer.time();

          std::cout << "Scales computed in " << time << " second(s) using "
                    << (memory>>20) << " MiB of memory:" << std::endl;
          std::cout << " * Global K scale: " << k_scale << std::endl;
          std::cout << " * Global range scale: " << range_scale << std::endl;
        }
      else
        {
          stream.seekg(0);
          std::cerr << "TRYING 2D CASE" << std::endl;
          std::vector<Point_2> points_2;

          // Read ASCII 2D point set file
          std::string str;
          while (getline (stream, str))
            {
              std::istringstream iss (str);
              double x = 0., y = 0.;
              iss >> x >> y;
              points_2.push_back (Point_2 (x, y));
            }

          if (points_2.empty())
            {
              std::cerr << "Error: cannot read file " << fname << std::endl;
              return EXIT_FAILURE;
            }
          
          task_timer.start();

          // estimate global scale
          std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale (points_2.begin(), points_2.end());
          FT range_scale = CGAL::estimate_global_range_scale (points_2.begin(), points_2.end());

          std::size_t memory = CGAL::Memory_sizer().virtual_size();
          double time = task_timer.time();

          std::cout << "Scales computed in " << time << " second(s) using "
                    << (memory>>20) << " MiB of memory:" << std::endl;
          std::cout << " * Global K scale: " << k_scale << std::endl;
          std::cout << " * Global range scale: " << range_scale << std::endl;
        }
    }
  else
    {
      std::cerr << "Error: cannot read file " << fname << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

