//
//file:  examples/Interpolation/sibson_interoplation_2.C 
//
#include <CGAL/basic.h>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>
#include <CGAL/interpolation_functions.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K>                Delaunay_triangulation;
typedef CGAL::Interpolation_gradient_fitting_traits_2<K> Traits;

typedef K::FT                                            Coord_type;
typedef K::Point_2                                       Point;
typedef std::map<Point, Coord_type, K::Less_xy_2>        Point_value_map ;
typedef std::map<Point, K::Vector_2 , K::Less_xy_2 >     Point_vector_map;

//Functor class for accessing the function values/gradients
template< class Map >  
struct DataAccess : public std::unary_function< typename Map::key_type,
		    typename Map::mapped_type> {
  typedef typename Map::mapped_type Data_type;
  typedef typename Map::key_type  Point;
  
  DataAccess< Map >(const Map& m): map(m){};
  
  Data_type operator()(const Point& p) { 
    typename Map::const_iterator mit = map.find(p);
    if(mit!= map.end())
      return mit->second;
    return Data_type();
  };
  
  const Map& map;
};

int main()
{
  
  Delaunay_triangulation T;

  Point_value_map values;
  Point_vector_map gradients;

  //parameters for spherical function:
  Coord_type a(0.25), bx(1.3), by(-0.7), c(0.2);
  for (int y=0 ; y<4 ; y++)
    for (int x=0 ; x<4 ; x++){ 
      K::Point_2 p(x,y);
      T.insert(p);
      values.insert(std::make_pair(p,a + bx* x+ by*y + c*(x*x+y*y)));
    }
  sibson_gradient_fitting_nn_2(T,std::inserter(gradients,gradients.begin()),
			       DataAccess<Point_value_map>(values), Traits());
  
 
  //coordiante computation
  K::Point_2 p(1.6,1.4);
  std::vector< std::pair< Point, Coord_type > > coords;
  Coord_type norm = 
    CGAL::natural_neighbor_coordinates_2(T, p,std::back_inserter(coords)).second;
  
  
  //Sibson interpolant: version without sqrt:
  Coord_type res =  CGAL::sibson_c1_interpolation_square(coords.begin(),
							 coords.end(),norm,p, 
							 DataAccess<Point_value_map>(values),
							 DataAccess<Point_vector_map>(gradients),
							 Traits());
  std::cout << "   Tested interpolation on " << p << " interpolation: " << res 
	    << " exact: " 
	    << a + bx * p.x()+ by * p.y()+ c*(p.x()*p.x()+p.y()*p.y()) 
	    << std::endl;
  return 1; 
}
