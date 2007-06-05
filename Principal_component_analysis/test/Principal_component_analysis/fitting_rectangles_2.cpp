// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

#include <vector>
#include <cassert>
#include <stdlib.h>
#define THRESHOLD 0.001
// types

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Iso_rectangle_2   Iso_rectangle_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Vector_2 Vector_2;

// case with only one square in container
// the fitting line must be y = 1/2
void test_1()
{
  
  std::list<Iso_rectangle_2> Iso_rectangles;
  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(0.0,0.0),Point_2(4.0,4.0)));
  //  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(4.0,0.0),Point_2(5.0,1.0)));

  // fit a line
  // call all versions of the function
  std::cout<<"Test 1"<<std::endl;
  std::cout << "fit 2D line...";
  Kernel k;
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,centroid,k,CGAL::PCA_dimension_2_tag());
  std::cout << "done (quality: " << quality << ") Line: " << line<<std::endl;

  
  if(!(line.is_horizontal() && std::abs(line.c()/line.b()- -2) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

// case with only one square vs two rectangles that it is composed of
// in container
// the fitting line must be the same
void test_2()
{
  std::cout<<"Test 2"<<std::endl;
  std::list<Iso_rectangle_2> Iso_rectangles;
  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(1.6,15.2),Point_2(5.6,19.2)));
  //  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(4.0,0.0),Point_2(5.0,1.0)));

  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line to bigger rectangle...";
  Kernel k;
  Line_2 line;
  Point_2 centroid;
  FT quality;
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,centroid,k,CGAL::PCA_dimension_2_tag());
  std::cout << "done (quality: " << quality << ") Line: " << line<<" centroid: "<<centroid<<std::endl;


  // fit a line
  // call all versions of the function
  std::list<Iso_rectangle_2> Iso_rectangles1;
  Iso_rectangles1.push_back(Iso_rectangle_2(Point_2(1.6,15.2),Point_2(3.6,19.2)));
  Iso_rectangles1.push_back(Iso_rectangle_2(Point_2(3.6,15.2),Point_2(5.6,19.2)));
  std::cout << "fit 2D line to two small rectangles...";
  Kernel k1;
  Line_2 line1;
  Point_2 centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_2(Iso_rectangles1.begin(),Iso_rectangles1.end(),line1,CGAL::PCA_dimension_2_tag());
  quality1 = linear_least_squares_fitting_2(Iso_rectangles1.begin(),Iso_rectangles1.end(),line1,centroid1,CGAL::PCA_dimension_2_tag());
  quality1 = linear_least_squares_fitting_2(Iso_rectangles1.begin(),Iso_rectangles1.end(),line1,centroid1,k1,CGAL::PCA_dimension_2_tag());
  std::cout << "done (quality: " << quality1 << ") Line: " << line1<<" centroid: "<<centroid1<<std::endl;

  
  if(!(std::abs(-1.0*line.a()/line.b() - -1.0*line1.a()/line1.b()) <= THRESHOLD && std::abs(line.c()/line.b() - line1.c()/line1.b()) <= THRESHOLD && std::abs(quality1 - quality) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

/*
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
  Kernel k;
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid);
  quality = linear_least_squares_fitting_2(points.begin(),points.end(),line,centroid,k);

  std::cout << "done (quality: " << quality << ")" <<" line: "<<line<< std::endl;

  std::list<Segment_2> segments;
  segments.push_back(Segment_2(p,q));


  // fit a line
  // call all versions of the function
  std::cout << "fit 2D line to segment...";
  Kernel k1;
  Line_2 line1;
  Point_2 centroid1;
  FT quality1;
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1);
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1,centroid1);
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1,centroid1,k1);
  quality1 = linear_least_squares_fitting_2(segments.begin(),segments.end(),line1,centroid1,k1,false);
  std::cout << "done (quality: " << quality1 << ")" <<" line: "<<line1<< std::endl;

  if(!(std::abs(-1.0*line.a()/line.b() - -1.0*line1.a()/line1.b()) <= THRESHOLD && std::abs(line.c()/line.b() - line1.c()/line1.b()) <= THRESHOLD && std::abs(quality1 - quality) <= THRESHOLD))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }

}

*/
int main()
{
  std::cout << "Test 2D linear_least_squares_fitting_rectangles"  << std::endl;

  test_1();
  test_2();


  return 0; // success
}
