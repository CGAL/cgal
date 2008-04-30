// outliers_removal_test.cpp


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the outliers removal methods
// No output
// Input files are .xyz
//----------------------------------------------------------
// outliers_removal_test points1.xyz points2.xyz...


// CGAL
#include <CGAL/basic.h> // include basic.h before testing #defines
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/remove_outliers_wrt_avg_knn_sq_distance_3.h>

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

// read point set from .xyz file
bool read_point_set(char *input_filename,
                    std::vector<Point>& points)
{
  std::cerr << "  Open " << input_filename << " for reading...";

  std::ifstream stream(input_filename);
  if(!stream.is_open())
  {
    std::cerr << "failed" << std::endl;
    return false;
  }

  // read point set
  Point point;
  while(!stream.fail())
  {
    stream >> point;
    points.push_back(point);
  }
  std::cerr << "ok (" << points.size() << " points)" << std::endl;
  return true;
}

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
    if(read_point_set(argv[i],points))
    {
      test_avg_knn_sq_distance(points, k, outliers_percentage);
    }
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}
