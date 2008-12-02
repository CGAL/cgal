// Example program for the linear_least_square_fitting function
// on a set of 3D tetrahedra
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
	Point a(1.0,2.0,3.0);
	Point b(4.0,5.0,6.0);
	Point c(7.0,8.0,9.0);
	Point d(8.0,7.0,6.0);
	Point e(2.0,4.0,5.0);

	std::list<Tetrahedron> tetrahedra;
  tetrahedra.push_back(Tetrahedron(a,b,c,d));
  tetrahedra.push_back(Tetrahedron(a,b,c,e));
  tetrahedra.push_back(Tetrahedron(b,c,d,e));

  Line line;
  Plane plane;

	// fit line and plane to whole tetrahedra
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<3>());

	// fit line and plane to tetrahedron faces
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<2>());
    
	// fit line and plane to tetrahedron edges
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<1>());

	// fit line and plane to tetrahedron vertices
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<0>());

	return 0;
}
