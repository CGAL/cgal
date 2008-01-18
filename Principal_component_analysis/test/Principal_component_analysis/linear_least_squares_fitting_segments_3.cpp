// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <CGAL/Testsuite/assert.h>
#include <cstdlib>
#define THRESHOLD 0.001
// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_3 Line_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Direction_3 Direction_3;

// case with only one segment in container
// the fitting line must be the same
void test_1()
{
  std::list<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(3.0,3.0,3.0),Point_3(1.0,1.0,1.0)));


  // fit a line
  // call all versions of the function
  std::cout<<"Test 1"<<std::endl;
  std::cout << "fit 3D line...";
  Kernel k;
  Line_3 line;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  std::cout << "done (quality: " << quality << ") line: " <<line<< std::endl;

  if(!(line.has_on(Point_3(2.0,2.0,2.0)) && (double)line.to_vector().x() - 1/std::sqrt(3.0) <=THRESHOLD && (double)line.to_vector().y() - 1.0/std::sqrt(3.0) <= THRESHOLD && (double)line.to_vector().z() - 1.0/std::sqrt(3.0) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


// case with two segments cutting the one segment in two parts in 
// container the fitting line must be the same
void test_2()
{
  std::cout<<"Test 2"<<std::endl;
  std::list<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(2.0,2.0,2.0),Point_3(1.5,1.5,1.5)));
  segments.push_back(Segment_3(Point_3(1.0,1.0,1.0),Point_3(1.5,1.5,1.5)));


  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line...";
  Kernel k;
  Line_3 line;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  std::cout << "done (quality: " << quality << ")" << std::endl;

  if(!(line.has_on(Point_3(2.0,2.0,2.0)) && (double)line.to_vector().x() - 1/std::sqrt(3.0) <=THRESHOLD && (double)line.to_vector().y() - 1.0/std::sqrt(3.0) <= THRESHOLD && (double)line.to_vector().z() - 1.0/std::sqrt(3.0) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

// case with four segments cutting the one segment in four parts in 
// container the fitting line must be the same
void test_3()
{
  std::cout<<"Test 3"<<std::endl;
  std::list<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(2.0,2.0,2.0),Point_3(1.75,1.75,1.75)));
  segments.push_back(Segment_3(Point_3(1.75,1.75,1.75),Point_3(1.50,1.50,1.50)));
  segments.push_back(Segment_3(Point_3(1.50,1.50,1.50),Point_3(1.25,1.25,1.25)));
  segments.push_back(Segment_3(Point_3(1.25,1.25,1.25),Point_3(1.0,1.0,1.0)));

  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line...";
  Kernel k;
  Line_3 line;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  std::cout.precision(20);
  std::cout << "done (quality: " << quality << ")"  <<" line: "<<line<<std::endl;

  if(!(line.point(0)==Point_3(1.5,1.5,1.5)  && (double)line.to_vector().x() - 1/std::sqrt(3.0) <=THRESHOLD && (double)line.to_vector().y() - 1.0/std::sqrt(3.0) <= THRESHOLD && (double)line.to_vector().z() - 1.0/std::sqrt(3.0) <= THRESHOLD))
  {
    std::cout << "failure " <<std::endl;
    exit(1); // failure Direction_3(0.57735,0.57735,0.57735)
  }
}


// case with a segments in container against just the two end points
// the fitting line must be the same
void test_4()
{
  std::cout<<"Test 4"<<std::endl;
  Point_3 p = Point_3(1.5,6.0,12.34);
  Point_3 q = Point_3(1.12,7.21,4.3);
  std::list<Segment_3> segments;
  segments.push_back(Segment_3(p,q));


  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line to segment...";
  Kernel k;
  Line_3 line;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  std::cout << "done (quality: " << quality << ")" <<" line: "<<line<< std::endl;

  std::list<Point_3> points;
  points.push_back(p);
  points.push_back(q);

  // fit a line
  // call all versions of the function
  std::cout << "fit 3D line to end points...";
  Kernel k1;
  Line_3 line1;
  Point_3 centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_3(points.begin(),points.end(),line1,CGAL::PCA_dimension_0_tag());
  quality1 = linear_least_squares_fitting_3(points.begin(),points.end(),line1,centroid1,CGAL::PCA_dimension_0_tag());
  std::cout << "done (quality: " << quality1 << ")" <<" line: "<<line1<<std::endl;

  if(!(line.has_on(line1.point(0)) && (double)line.to_vector().y()/line.to_vector().x() - (double)line1.to_vector().y()/line1.to_vector().x() <= THRESHOLD && (double)line.to_vector().z()/line.to_vector().x() - (double)line1.to_vector().z()/line1.to_vector().x() <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}


//:::::PLANE FITTING:::::::::::
// case with only one segment in container
// the fitting plane must be the same
void test_5()
{
  std::vector<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(3.0,3.0,3.0),Point_3(1.0,1.0,1.0)));


  // fit a plane
  // call all versions of the function
  std::cout<<"Test 5"<<std::endl;
  std::cout << "fit 3D plane...";
  Kernel k;
  Plane_3 plane;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::PCA_dimension_1_tag());
  std::cout << "done (quality: " << quality << ") plane: " <<plane<< std::endl;

  if(!plane.has_on(segments[0].supporting_line()))
  //  if(!(line.has_on(Point_3(2.0,2.0,2.0)) && (double)line.to_vector().x() - 1/std::sqrt(3.0) <=THRESHOLD && (double)line.to_vector().y() - 1.0/std::sqrt(3.0) <= THRESHOLD && (double)line.to_vector().z() - 1.0/std::sqrt(3.0) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

void test_6() {
  std::vector<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(3.0,3.0,3.0),Point_3(1.0,1.0,1.0)));


  // call all versions of the function
  std::cout<<"Test 6"<<std::endl;
  Kernel k;
  Plane_3 plane;
  Line_3 line;
  Point_3 centroid;
  FT quality;
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::PCA_dimension_0_tag());

  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,centroid,CGAL::PCA_dimension_0_tag());

  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,CGAL::PCA_dimension_0_tag());

  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  quality = linear_least_squares_fitting_3(segments.begin(),segments.end(),line,centroid,CGAL::PCA_dimension_0_tag());
}
int main()
{
  std::cout << "Test linear least squares fitting of 3D segments"  << std::endl;
  test_1();
  test_2();
  test_3();  
  test_4();
  test_5();
  test_6();
  return 0; // success
}
