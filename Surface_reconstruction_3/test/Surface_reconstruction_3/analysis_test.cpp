// analysis_test.cpp

// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the analysis methods
// Input files are .xyz
// No output
//----------------------------------------------------------
// analysis_test points1.xyz points2.xyz...


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/average_spacing_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

// STL stuff
#include <deque>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

void test_average_spacing(std::deque<Point>& points, // input point set
                          unsigned int nb_neighbors_avg_spacing) // number of neighbors
{
  std::cerr << "Compute average spacing to knn (knn="<< nb_neighbors_avg_spacing << ")... ";
  CGAL::Timer task_timer; task_timer.start();

  typedef std::deque<Point>::iterator Iterator;
  std::cerr << CGAL::average_spacing_3<Iterator,FT>(points.begin(), points.end(),
                                                    nb_neighbors_avg_spacing) << "\n";
  
  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "ok: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Analysis test" << std::endl;
  std::cerr << "No output" << std::endl;

  // decode parameters
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
    return EXIT_FAILURE;
  }

  // Average Spacing options
  const unsigned int nb_neighbors_avg_spacing = 7; // K-nearest neighbors = 1 ring (average spacing)
                                                   // LS: was 8

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for(int i=1; i<argc; i++)
  {
    std::cerr << std::endl;

    // Load point set
    std::deque<Point> points;
    std::cerr << "Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(argv[i], 
                                             std::back_inserter(points), 
                                             false /*skip normals*/))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
      
      test_average_spacing(points,nb_neighbors_avg_spacing);
    }
    else
    {
      std::cerr << "Error: cannot read file " << argv[i] << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
    }
  } // for each input file

  std::cerr << std::endl;

  // Return accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}
 
