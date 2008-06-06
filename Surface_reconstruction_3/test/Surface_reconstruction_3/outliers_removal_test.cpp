// outliers_removal_test.cpp

// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the outliers removal methods
// Input files are .xyz
// No output
//----------------------------------------------------------
// outliers_removal_test points1.xyz points2.xyz...


#include <CGAL/basic.h> // include basic.h before testing #defines

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/remove_outliers_wrt_avg_knn_sq_distance_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

// STL stuff
#include <list>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iterator>

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

void test_avg_knn_sq_distance(const std::vector<Point>& points, // input point set
                              unsigned int k,// number of neighbors
                              double outliers_percentage) // percentage of points to remove 
{
  std::cerr << "  Remove outliers using Average KNN Squared Distance...";

  // todo: use mutating version when ready
  std::vector<Point> output;
  CGAL::remove_outliers_wrt_avg_knn_sq_distance_3(
          points.begin(), points.end(),
          std::back_inserter(output),
          k, 
          outliers_percentage);
  
  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Normal estimation test" << std::endl;
  std::cerr << "No output" << std::endl;

  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz ..." << std::endl;
    return EXIT_FAILURE;
  }

  // for each input file
  const unsigned int k = 10; // # neighbors
  const unsigned int outliers_percentage = 5; // percentage of points to remove
  for(int i=1; i<argc; i++)
  {
    std::vector<Point> points;
    std::cerr << "  Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(argv[i], 
                                             std::back_inserter(points), 
                                             false /*skip normals*/))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
      
      test_avg_knn_sq_distance(points, k, outliers_percentage);
    }
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}
 
