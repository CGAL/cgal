// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/algorithm.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

#include <vector>
#include <cassert>
#include <cstdlib>
#define THRESHOLD 0.001
// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Vector_2 Vector_2;

// case with only one segment in container
// the fitting line must be the same
void test_1()
{
  std::list<Segment_2> segments;
  segments.push_back(Segment_2(Point_2(2.0,2.0),Point_2(1.0,1.0)));


  // fit a line
  // call all versions of the function
  std::cout<<"Test 1"<<std::endl;
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - 1) <= THRESHOLD && std::abs(line.c()/line.b()) <= THRESHOLD && 1 - quality <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}

// case with two segments cutting the one segment in two parts in 
// container the fitting line must be the same
void test_2()
{
  std::cout<<"Test 2"<<std::endl;
  std::list<Segment_2> segments;
  segments.push_back(Segment_2(Point_2(2.0,2.0),Point_2(1.5,1.5)));
  segments.push_back(Segment_2(Point_2(1.0,1.0),Point_2(1.5,1.5)));


  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - 1) <= THRESHOLD && std::abs(line.c()/line.b()) <= THRESHOLD && 1 - quality <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}

// case with four segments cutting the one segment in four parts in 
// container the fitting line must be the same
void test_3()
{
  std::cout<<"Test 3"<<std::endl;
  std::list<Segment_2> segments;
  segments.push_back(Segment_2(Point_2(2.0,2.0),Point_2(1.75,1.75)));
  segments.push_back(Segment_2(Point_2(1.0,1.0),Point_2(1.25,1.25)));
  segments.push_back(Segment_2(Point_2(1.5,1.5),Point_2(1.25,1.25)));
  segments.push_back(Segment_2(Point_2(1.5,1.5),Point_2(1.75,1.75)));


  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line...";
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - 1) <= THRESHOLD && std::abs(line.c()/line.b()) <= THRESHOLD && 1 - quality <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}

// case with a segments in container against just the two end points
// the fitting line must be the same
void test_4()
{
  std::cout<<"Test 4"<<std::endl;
  Point_2 p = Point_2(1.0,0.0);
  Point_2 q = Point_2(0.0,1.0);
  std::list<Segment_2> segments;
  segments.push_back(Segment_2(p,q));


  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line to segment...";
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  quality = linear_least_squares_fitting_2(segments.begin(),segments.end(),line,centroid,CGAL::Dimension_tag<1>());
  std::cout << "done (quality: " << quality << ")" <<" line: "<<line<< std::endl;

  std::list<Point_2> points;
  points.push_back(p);
  points.push_back(q);

  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line to end points...";
  Line_2 line1;
  Point_2 centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_2(points.begin(),points.end(),line1,CGAL::Dimension_tag<0>());
  quality1 = linear_least_squares_fitting_2(points.begin(),points.end(),line1,centroid1,CGAL::Dimension_tag<0>());
  std::cout << "done (quality: " << quality1 << ")" <<" line: "<<line1<<std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - -1.0*line1.a()/line1.b()) <= THRESHOLD && std::abs(line.c()/line.b() - line1.c()/line1.b()) <= THRESHOLD && std::abs(quality1 - quality) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }
}


// case with a point set on a segment vs the segment
// the fitting line must be the same segment
void test_5(const unsigned int nb_points)
{
  std::cout<<"Test 5"<<std::endl;
  // create points on a horizontal segment
  Point_2 p(2.0,0.0);
  Point_2 q(5.0,5.0);

  std::cout << "generate " << nb_points << 
       " 2D points on a horizontal line...";
  std::list<Point_2> points;
  points_on_segment_2(p,q,100,std::back_inserter(points));
  std::cout << "done " << std::endl;

  // fit a line
  std::cout << "fit 2D line to points...";
  Line_2 line;
  Point_2 centroid;

  // call all versions of the function
  FT quality;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,CGAL::Dimension_tag<0>());

  std::cout << "done (quality: " << quality << ")" <<" line: "<<line<< std::endl;

  std::list<Segment_2> segments;
  segments.push_back(Segment_2(p,q));


  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line to segment...";
  Line_2 line1;
  Point_2 centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1,CGAL::Dimension_tag<1>());
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1,centroid1,CGAL::Dimension_tag<1>());
  std::cout << "done (quality: " << quality1 << ")" <<" line: "<<line1<< std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - -1.0*line1.a()/line1.b()) <= THRESHOLD && std::abs(line.c()/line.b() - line1.c()/line1.b()) <= THRESHOLD && std::abs(quality1 - quality) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    std::exit(1); // failure
  }

}


int main()
{
  std::cout << "Test 2D linear least squares fitting of segments"  << std::endl;

  test_1();
  test_2();
  test_3();  
  test_4();
  test_5(20);
  //  test_2D_point_set(100);

  return 0; // success
}
