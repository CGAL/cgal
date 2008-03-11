// Example program for the linear_least_square_fitting function
// on a set of tetrahedra in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>
#include <cstdlib> // for std::rand

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Tetrahedron_3     Tetrahedron;

FT random_value()
{
	return (FT)std::rand() / (FT)RAND_MAX;
}

Point random_point()
{
	return Point(random_value(),
		           random_value(),
							 random_value());
}

int main(void)
{
	Point a = random_point();
	Point b = random_point();
	Point c = random_point();
	Point d = random_point();
	Point e = random_point();

	std::list<Tetrahedron> tetrahedra;
  tetrahedra.push_back(Tetrahedron(a,b,c,d));
  tetrahedra.push_back(Tetrahedron(a,b,c,e));

  Line line;
  Plane plane;

	// fit a line and a plane to tetrahedra
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::PCA_dimension_3_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_3_tag());

	// fit a line and a plane to tetrahedron faces
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::PCA_dimension_2_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_2_tag());
    
	// fit a line and a plane to tetrahedron edges
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::PCA_dimension_1_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_1_tag());

	// fit a line and a plane to tetrahedron vertices
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line, CGAL::PCA_dimension_0_tag());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::PCA_dimension_0_tag());

	return 0;
}
