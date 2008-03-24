#include <CGAL/basic.h> // include basic.h before testing #defines

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>

#include <CGAL/estimate_normals_pca_3.h>
#include <CGAL/estimate_normals_jet_fitting_3.h>
#include <CGAL/Oriented_normal_3.h>

// STL stuff
#include <list>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iterator>

// types
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Oriented_normal_3<Kernel> Normal;

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

void test_pca(std::list<Point>& points,
							const unsigned int k)
{
	std::cerr << "  Estimate normals using KNN and point-based PCA...";
	std::list<Normal> normals;
	CGAL::estimate_normals_pca_3(points.begin(),points.end(),std::back_inserter(normals),k);
	std::cerr << "ok" << std::endl;
}

void test_jet_fitting(std::list<Point>& points,
							        const unsigned int k)
{
	std::cerr << "  Estimate normals using KNN and jet fitting...";
	std::list<Normal> normals;
	CGAL::estimate_normals_jet_fitting_3(points.begin(),points.end(),std::back_inserter(normals),k);
	std::cerr << "ok" << std::endl;
}

// main function
int main(int argc, 
				 char * argv[])
{
  std::cerr << "Normal estimation test" << std::endl;

  if(argc < 2)
  {
		std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz ..." << std::endl;
		return EXIT_FAILURE;
  }

  // for each input file
	const unsigned int k = 10; // # neighbors
  for(int i=1; i<argc; i++)
  {
		std::list<Point> points;
		if(read_point_set(argv[i],points))
		{
			test_pca(points,k);
			test_jet_fitting(points,k);
		}
		else
			std::cerr << "  Unable to open file " << argv[i] << std::endl;
	}
	return EXIT_SUCCESS;
}



