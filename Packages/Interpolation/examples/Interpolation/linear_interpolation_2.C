//
//file: examples/Interpolation/linear_interoplation_2.C 
//
#include <CGAL/basic.h>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>               Traits;
typedef K::FT                                         Coord_type;

//Functor for accessing the function values
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
  std::map<Point, Coord_type, K::Less_xy_2> values;
  typedef DataAccess< std::map<Point, Coord_type, K::Less_xy_2 > > 
                                                        Value_access;
  
  Coord_type a(0.25), bx(1.3), by(-0.7);

  for (int y=0 ; y<3 ; y++)
    for (int x=0 ; x<3 ; x++){ 
      K::Point_2 p(x,y);
      T.insert(p);
      values.insert(std::make_pair(p,a + bx* x+ by*y));
    }
  //coordiante computation
  K::Point_2 p(1.3,0.34);
  std::vector< std::pair< Point, Coord_type > > coords;
  Coord_type norm = 
    CGAL::natural_neighbor_coordinates_2(T, p,std::back_inserter(coords)).second;  
 
  Coord_type res =  CGAL::linear_interpolation(coords.begin(), coords.end(), 
				               norm,Value_access(values));
  
  std::cout << "   Tested interpolation on " << p << " interpolation: " << res 
	    << " exact: " << a + bx* p.x()+ by* p.y()<< std::endl;
  return 1; 
}
