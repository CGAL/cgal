// Example program for the linear_least_square_fitting function
// on a set of tetrahedra in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  Kernel;
typedef Kernel::Line_3            Line;
typedef Kernel::Plane_3           Plane;
typedef Kernel::Point_3           Point;
typedef Kernel::Tetrahedron_3     Tetrahedron;

int main()
{
	// generate two tetrahedra
	std::list<Tetrahedron> tetrahedra;
	Point a(0.0,0.0,0.0);
	Point b(1.0,0.0,0.0);
	Point c(0.0,1.0,0.0);
	Point d(0.0,0.0,1.0);
	Point e(0.0,1.0,1.0);
  tetrahedra.push_back(Tetrahedron(a,b,c,d));
  tetrahedra.push_back(Tetrahedron(a,b,c,e));

  Kernel kernel;
  Point centroid;
  Plane plane;
  Line line;

	// fit line, no centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_0_tag());

	// fit line, centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::PCA_dimension_0_tag());

	// fit plane, no centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_0_tag());

	// fit plane, centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::PCA_dimension_0_tag());

  return 0;
}
