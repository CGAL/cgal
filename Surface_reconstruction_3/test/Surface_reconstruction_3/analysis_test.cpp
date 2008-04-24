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
#include <CGAL/average_spacing_3.h>

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
                    std::list<Point>& points)
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
    if(read_point_set(argv[i],points))
    {
      test_average_spacing(points,k);
    }
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}
 

