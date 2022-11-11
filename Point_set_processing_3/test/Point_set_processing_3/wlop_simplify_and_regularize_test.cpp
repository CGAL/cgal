// wlop_simplify_and_regularize_test.cpp

//----------------------------------------------------------
// Test the wlop simplify and regularize method:
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// wlop_simplify_and_regularize_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/tags.h>

// This package
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <string>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Removes outliers
template<typename Concurrency_tag>
void test_wlop_simplify_and_regularize(
                              std::vector<Point>& points, // input point set
                              std::vector<Point>& output,
                              double retain_percentage, // percentage of points to remove
                              double neighbor_radius, // neighborhood size
                              unsigned int iter_number, // iteration number
                              bool need_compute_density)

{
  CGAL::Real_timer task_timer; task_timer.start();
  std::cerr << "Running WLOP simplify and regularize, (retain_percentage: "
            << retain_percentage << "%, neighbor_radius="
            << neighbor_radius << ")...\n";

  // Make room for sample points
  std::vector<Point> points_sampled;
  points_sampled.resize(static_cast<std::size_t>(points.size() * (retain_percentage / 100.)));

  output.clear();
  // Run algorithm
  CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>
    (points, std::back_inserter(output),
     CGAL::parameters::select_percentage(retain_percentage).
     neighbor_radius(neighbor_radius).
     number_of_iterations(iter_number).
     require_uniform_sampling(need_compute_density));

  output.clear();

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "ok: " << task_timer.time() << " seconds, "
                      << (memory>>20) << " Mb allocated"
                      << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "WLOP simplify and regularize" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "For each input point set, apply WLOP algorithm.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  //Algorithm parameters
  const double retain_percentage = 2;   // percentage of points to retain.
  const double neighbor_radius = 0.5;   // neighbors size.
  const unsigned int iter_number = 25;     // number of iterations.
  const bool need_compute_density = false;  // if needed to compute density.

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for(int i=1; i<argc; i++)
  {
    std::cerr << std::endl;

    //***************************************
    // Loads point set
    //***************************************

    // Reads the point set file in points[].
    std::vector<Point> points;
    std::cerr << "Opening " << argv[i] << " for reading..." << std::endl;

    // If XYZ file format:
    if(CGAL::IO::read_points(argv[i], std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << argv[i] << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    //***************************************
    // Test
    //***************************************
#ifdef CGAL_LINKED_WITH_TBB
    std::vector<Point> points2(points);
#endif
    std::vector<Point> output;
    test_wlop_simplify_and_regularize<CGAL::Sequential_tag>(
      points, output, retain_percentage, neighbor_radius,
      iter_number, need_compute_density);

#ifdef CGAL_LINKED_WITH_TBB
    output.clear();
    test_wlop_simplify_and_regularize<CGAL::Parallel_tag>(
      points2, output, retain_percentage, neighbor_radius,
      iter_number, need_compute_density);
#endif // CGAL_LINKED_WITH_TBB


  } // for each input file

  std::cerr << std::endl;
  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}
