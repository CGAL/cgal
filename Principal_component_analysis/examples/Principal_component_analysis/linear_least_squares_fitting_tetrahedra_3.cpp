// Example program for the linear_least_square_fitting function
// on a set of tetrahedra in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3           Plane_3;
typedef K::Point_3           Point_3;
typedef K::Tetrahedron_3     Tetrahedron_3;

int main(void)
{
	Point a( 0.0,0.0,0.0);
	Point b( 1.0,0.0,0.0);
	Point c(-1.0,0.0,0.0);
	Point d( 0.0,1.0,1.0);
	Point e( 0.0,0.0,0.0);

	std::list<Tetrahedron_3> tetrahedra;
  tetrahedra.push_back(Tetrahedron_3(a,b,c,d));
  tetrahedra.push_back(Tetrahedron_3(a,b,c,e));

	// fit a line
  Line_3 line;
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::PCA_dimension_3_tag());

	// fit a plane
  Plane_3 plane;
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_3_tag());
    
  return 0;
}
