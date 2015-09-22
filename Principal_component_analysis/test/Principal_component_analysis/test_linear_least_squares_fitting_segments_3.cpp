// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/algorithm.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Default_diagonalize_traits.h>

#include <vector>
#include <cassert>
#include <cstdlib>
#define THRESHOLD 0.001
// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_3 Line;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Point_3 Point;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Direction_3 Direction;

FT fit_set(std::list<Segment>& segments,
           Plane& plane,
           Line& line)
{
  // fit a plane
  // call all versions of the function
  Kernel kernel;
  FT quality;
  Point centroid;

  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>(),kernel,
					   CGAL::Default_diagonalize_traits<FT,3>());

  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::Dimension_tag<1>(),kernel,
					   CGAL::Default_diagonalize_traits<FT,3>());

  return quality;
}


// case with only one segment in container
// the fitting line must be the same
void test_1()
{
  std::list<Segment> segments;
  segments.push_back(Segment(Point(3.0,3.0,3.0),Point(1.0,1.0,1.0)));

  // fit a line
  // call all versions of the function
  std::cout << "Test 1" << std::endl;
  std::cout << "fit 3D line...";
  Line line;
	Plane plane;
  Point centroid;

	FT quality = fit_set(segments,plane,line);
  std::cout << "done, quality: " << quality << std::endl;

  if(!line.has_on(Point(3.0,3.0,3.0)) || 
		 !line.has_on(Point(1.0,1.0,1.0)))
  {
    std::cout << "failure" << std::endl;
		std::cout << "line: " << line << std::endl;
    // exit(1); // failure
  }
}


// case with two segments cutting the one segment in two parts in 
// container the fitting line must be the same
void test_2()
{
  std::cout<<"Test 2"<<std::endl;
  std::list<Segment> segments;
  segments.push_back(Segment(Point(1.0,1.0,1.0),Point(1.5,1.5,1.5)));
  segments.push_back(Segment(Point(2.0,2.0,2.0),Point(1.5,1.5,1.5)));

  // fit a line
  std::cout << "fit 3D line...";
  Line line;
	Plane plane;
  Point centroid;
	FT quality = fit_set(segments,plane,line);
  std::cout << "done, quality: " << quality << std::endl;

  if(!line.has_on(Point(2.0,2.0,2.0)) ||
		 !line.has_on(Point(1.0,1.0,1.0)))
  {
    std::cout << "failure" << std::endl;
		std::cout << "line: " << line << std::endl;
    // exit(1); // failure
  }
}

// case with four segments cutting the one segment in four parts in 
// container the fitting line must be the same
void test_3()
{
  std::cout<<"Test 3"<<std::endl;
  std::list<Segment> segments;
  segments.push_back(Segment(Point(2.00,2.00,2.00),Point(1.75,1.75,1.75)));
  segments.push_back(Segment(Point(1.75,1.75,1.75),Point(1.50,1.50,1.50)));
  segments.push_back(Segment(Point(1.50,1.50,1.50),Point(1.25,1.25,1.25)));
  segments.push_back(Segment(Point(1.25,1.25,1.25),Point(1.00,1.00,1.00)));

  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line...";
  Line line;
	Plane plane;
  Point centroid;
	FT quality = fit_set(segments,plane,line);
  std::cout.precision(20);
  std::cout << "done, quality: " << quality << std::endl;

  if(!line.has_on(Point(2.0,2.0,2.0)) ||
		 !line.has_on(Point(1.0,1.0,1.0)))
  {
    std::cout << "failure" << std::endl;
		std::cout << "line: " << line << std::endl;
    // exit(1); // failure
  }
}


// case with a segments in container against just the two end points
// the fitting line must be the same
void test_4()
{
  std::cout<<"Test 4"<<std::endl;
  Point p = Point(1.5,6.0,12.34);
  Point q = Point(1.12,7.21,4.3);
  std::list<Segment> segments;
  segments.push_back(Segment(p,q));


  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line to segment...";
  Line line;
	Plane plane;
  Point centroid;
	FT quality = fit_set(segments,plane,line);
  std::cout << "done, quality: " << quality << std::endl;

  std::list<Point> points;
  points.push_back(p);
  points.push_back(q);

  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line to end points...";
  Line line1;
  Point centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_3(points.begin(),points.end(),line1,CGAL::Dimension_tag<0>());
  quality1 = linear_least_squares_fitting_3(points.begin(),points.end(),line1,centroid1,CGAL::Dimension_tag<0>());
  std::cout << "done, quality: " << quality1 << std::endl;

  if(!(line.has_on(line1.point(0)) && (double)line.to_vector().y()/line.to_vector().x() -
		(double)line1.to_vector().y()/line1.to_vector().x() <= THRESHOLD &&
		(double)line.to_vector().z()/line.to_vector().x() -
		(double)line1.to_vector().z()/line1.to_vector().x() <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    // exit(1); // failure
  }
}


// Plane fitting
// case with only one segment in container
// the fitting plane must be the same
void test_5()
{
  std::vector<Segment> segments;
  segments.push_back(Segment(Point(3.0,3.0,3.0),Point(1.0,1.0,1.0)));


  // fit a plane
  // call all versions of the function
  std::cout<<"Test 5"<<std::endl;
  std::cout << "fit 3D plane...";
	Line line;
	Plane plane;
  Point centroid;
  FT quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  std::cout << "done, quality: " << quality << std::endl;

  if(!plane.has_on(segments[0].supporting_line()))
  {
    std::cout << "failure" << std::endl;
    // exit(1); // failure
  }
}

void test_6()
{
  std::vector<Segment> segments;
  segments.push_back(Segment(Point(3.0,3.0,3.0),Point(1.0,1.0,1.0)));

  // call all versions of the function
  std::cout << "Test 6" << std::endl;
  Plane plane;
  Line line;
  Point centroid;

	// fit plane
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::Dimension_tag<0>());

	// fit line
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<0>());
}

int main()
{
  std::cout << "Test 3D linear least squares fitting of segments"  << std::endl;
  test_1();
  test_2();
  test_3();  
  test_4();
  test_5();
  test_6();
  return 0; // success
}
