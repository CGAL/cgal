// smoothing_test.cpp

// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the smoothing methods
// Input files are .xyz
// No output
//----------------------------------------------------------
// smoothing_test points1.xyz points2.xyz...


#include <CGAL/basic.h> // include basic.h before testing #defines

#ifdef CGAL_USE_LAPACK


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

#include <CGAL/smooth_jet_fitting_3.h>

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
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

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

void test_jet_fitting(std::list<Point>& points,
                      const unsigned int k)
{
  std::cerr << "  Smooth using KNN and jet fitting...";
  std::list<Point> output;
  CGAL::smooth_jet_fitting_3(points.begin(),points.end(),std::back_inserter(output),k);
  std::cerr << "ok" << std::endl;

  // mutating version of the same function
  std::cerr << "  Smooth using KNN and jet fitting...";
  CGAL::smooth_jet_fitting_3(points.begin(),points.end(),k);
  std::cerr << "ok" << std::endl;
}

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Smoothing test" << std::endl;

  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz ..." << std::endl;
    return EXIT_FAILURE;
  }

  // for each input file
  const unsigned int k = 20; // # neighbors
  for(int i=1; i<argc; i++)
  {
    std::list<Point> points;
    if(read_point_set(argv[i],points))
      test_jet_fitting(points,k);
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}
 

#else // CGAL_USE_LAPACK


#include <iostream>
#include <cstdlib>

// ----------------------------------------------------------------------------
// Empty main() if LAPACK is not installed
// ----------------------------------------------------------------------------

int main()
{
    std::cerr << "Skip smoothing test as LAPACK is not installed" << std::endl;
    return EXIT_SUCCESS;
}

#endif // CGAL_USE_LAPACK

