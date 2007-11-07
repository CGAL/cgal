// test for the linear_least_square_fitting() functions.


#include <CGAL/Cartesian.h>
#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <vector>
#include <CGAL/Testsuite/assert.h>
#include <stdlib.h>

// types
typedef CGAL::Cartesian<float> Kernel;
typedef Kernel::FT FT;

typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;

void fit_point_set(std::list<Point>& points,
                   Plane& plane,
                   Line& line)
{
  // fit a plane
  // call all versions of the function
  std::cout << "fit 3D plane...";
  Kernel k;
  FT quality;
  Point centroid;

  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,centroid,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),line,centroid,k,CGAL::PCA_dimension_0_tag());

  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid,CGAL::PCA_dimension_0_tag());
  quality = linear_least_squares_fitting_3(points.begin(),points.end(),plane,centroid,k,CGAL::PCA_dimension_0_tag());

  std::cout << "done (quality: " << quality << ")" << std::endl;
}

// case with only one point in container
// the fitting plane must be horizontal by default
void test_3D()
{
  std::list<Point> points;
  points.push_back(Point(0,0,0));

  // fit plane
  Plane plane;
  Line line;
  fit_point_set(points,plane,line);

  Plane horizontal_plane(Point(0,0,0),Vector(0,0,1));
  if(!parallel(horizontal_plane,plane))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

Point random_point_xy()
{
  FT x = (FT)((double)rand() / (double)RAND_MAX);
  FT y = (FT)((double)rand() / (double)RAND_MAX);
  return Point(x,y,0);
}

// case with a random point set on a horizontal plane
// the fitting plane must be horizontal
void test_3D_point_set(const unsigned int nb_points)
{
  std::list<Point> points;
  unsigned int i;
  for(i=0;i<nb_points;i++)
    points.push_back(random_point_xy());

  // fit plane
  Plane plane;
  Line line;
  fit_point_set(points,plane,line);

  Plane horizontal_plane(Point(0,0,0),Vector(0,0,1));
  if(!parallel(horizontal_plane,plane))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

// case with a point set on a horizontal plane
// the fitting plane must be horizontal
void test_3D_triangle_set(const unsigned int nb_triangles)
{
  std::list<Triangle> triangles;
  unsigned int i;
  for(i=0;i<nb_triangles;i++)
    triangles.push_back(Triangle(random_point_xy(),random_point_xy(),random_point_xy()));

  // fit a plane
  // call all versions of the function
  std::cout << "fit 3D plane...";
  Kernel k;
  Line line;
  FT quality;
  Plane plane;
  Point centroid;

  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line,centroid,k,CGAL::PCA_dimension_2_tag());

  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,CGAL::PCA_dimension_2_tag());
  quality = linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,centroid,k,CGAL::PCA_dimension_2_tag());

  std::cout << "done (quality: " << quality << ")" << std::endl;

  Plane horizontal_plane(Point(0,0,0),Vector(0,0,1));
  if(!parallel(horizontal_plane,plane))
  {
    std::cout << "failure" << std::endl;
    exit(1); // failure
  }
}

int main()
{
  std::cout << "Test 3D linear_least_squares_fitting"  << std::endl;

  // 3D
  test_3D();
  test_3D_point_set(100);
  test_3D_triangle_set(100);

  return 0; // success
}
