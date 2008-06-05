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


#include <CGAL/basic.h> // include basic.h before testing #defines

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

// This package
#include <CGAL/average_spacing_3.h>
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

void test_average_spacing(std::list<Point>& points,
                          const unsigned int k)
{
  std::cerr << "  Average spacing using KNN: ";
  std::list<Point> output;
  typedef std::list<Point>::iterator Iterator;
  std::cerr << CGAL::average_spacing_3<Iterator,FT>(points.begin(),points.end(),k);
  std::cerr << " ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Analysis tests" << std::endl;
  std::cerr << "No output" << std::endl;

  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz ..." << std::endl;
    return EXIT_FAILURE;
  }

  // for each input file
  const unsigned int k = 8; // # neighbors
  for(int i=1; i<argc; i++)
  {
    std::list<Point> points;
    std::cerr << "  Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(argv[i], std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
      
      test_average_spacing(points,k);
    }
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}
 
