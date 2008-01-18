// Example program for the linear_least_square_fitting function
// on a set of tetrahedra in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Tetrahedron_3     Tetrahedron;

int main(void)
{
	Point a( 0.0,0.0,0.0);
	Point b( 1.0,0.0,0.0);
	Point c(-1.0,0.0,0.0);
	Point d( 0.0,1.0,1.0);
	Point e( 0.0,0.0,0.0);

	std::list<Tetrahedron> tetrahedra;
  tetrahedra.push_back(Tetrahedron(a,b,c,d));
  tetrahedra.push_back(Tetrahedron(a,b,c,e));

	// fit a line
  Line line;
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_3_tag());

	// fit a plane
  Plane plane;
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_3_tag());
    
  return 0;
}
