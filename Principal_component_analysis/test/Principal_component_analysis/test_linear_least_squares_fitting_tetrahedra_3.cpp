// Example program for the linear_least_square_fitting function
// on a set of tetrahedra in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Diagonalize_traits.h>
#endif

typedef double                FT;
typedef CGAL::Cartesian<FT>   Kernel;
typedef Kernel::Line_3        Line;
typedef Kernel::Plane_3       Plane;
typedef Kernel::Point_3       Point;
typedef Kernel::Tetrahedron_3 Tetrahedron;

int main()
{
  std::cout << "Test 3D linear least squares fitting of tetrahedra"  << std::endl;

	// generate two tetrahedra
	std::list<Tetrahedron> tetrahedra;
	Point a(0.0,0.0,0.0);
	Point b(1.0,0.0,0.0);
	Point c(0.0,1.0,0.0);
	Point d(0.0,0.0,1.0);
	Point e(0.0,1.0,1.0);
  tetrahedra.push_back(Tetrahedron(a,b,c,d));
  tetrahedra.push_back(Tetrahedron(a,b,c,e));

  Line line;
  Plane plane;
  Point centroid;

	// fit line, no centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,CGAL::Dimension_tag<0>());

	// fit line, centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<0>());

	// fit plane, no centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,CGAL::Dimension_tag<0>());

	// fit plane, centroid
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<0>());

  // If Eigen is available, it's used by default everywhere.
  // These additional lines test the fallback version
#ifdef CGAL_EIGEN3_ENABLED
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<3>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),line,centroid,CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());

  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<3>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<2>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<1>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
  linear_least_squares_fitting_3(tetrahedra.begin(),tetrahedra.end(),plane,centroid,CGAL::Dimension_tag<0>(),
                                 Kernel(), CGAL::Diagonalize_traits<FT, 3>());
#endif

  return 0;
}
