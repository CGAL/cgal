#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>               Delaunay_triangulation;
typedef CGAL::Interpolation_gradient_fitting_traits_2<K> Traits;

typedef K::FT                                            Coord_type;
typedef K::Point_2                                       Point;
typedef std::map<Point, Coord_type, K::Less_xy_2>        Point_value_map ;
typedef std::map<Point, K::Vector_2 , K::Less_xy_2 >     Point_vector_map;

int main()
{
  Delaunay_triangulation T;

  Point_value_map function_values;
  Point_vector_map function_gradients;

  //parameters for spherical function:
  Coord_type a(0.25), bx(1.3), by(-0.7), c(0.2);
  for (int y=0 ; y<4 ; y++)
    for (int x=0 ; x<4 ; x++){
      K::Point_2 p(x,y);
      T.insert(p);
      function_values.insert(std::make_pair(p,a + bx* x+ by*y + c*(x*x+y*y)));
    }
  sibson_gradient_fitting_nn_2(T,std::inserter(function_gradients,
					       function_gradients.begin()),
			       CGAL::Data_access<Point_value_map>
			       (function_values),
			       Traits());


  //coordiante computation
  K::Point_2 p(1.6,1.4);
  std::vector< std::pair< Point, Coord_type > > coords;
  Coord_type norm =
    CGAL::natural_neighbor_coordinates_2(T, p,std::back_inserter
					 (coords)).second;


  //Sibson interpolant: version without sqrt:
  std::pair<Coord_type, bool> res =
    CGAL::sibson_c1_interpolation_square
    (coords.begin(),
     coords.end(),norm,p,
     CGAL::Data_access<Point_value_map>(function_values),
     CGAL::Data_access<Point_vector_map>(function_gradients),
     Traits());
  if(res.second)
    std::cout << "   Tested interpolation on " << p
	      << " interpolation: " << res.first << " exact: "
	      << a + bx * p.x()+ by * p.y()+ c*(p.x()*p.x()+p.y()*p.y())
	      << std::endl;
  else
    std::cout << "C^1 Interpolation not successful." << std::endl
	      << " not all function_gradients are provided."  << std::endl
	      << " You may resort to linear interpolation." << std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
