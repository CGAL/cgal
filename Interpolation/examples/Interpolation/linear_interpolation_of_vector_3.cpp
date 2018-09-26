#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>                   Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>                     Traits;
typedef K::Vector_3                                         Vector_3;
typedef K::Point_2                                          Point_2;


namespace boost {
  
template<class K>
class value_initialized<CGAL::Vector_3<K> >
{
  private :
  
    typedef CGAL::Vector_3<K> T;
    T m_data;

  public :
    
    value_initialized()
    :
      m_data(0,0,0)
    { }

    T const & data() const
    {
      return m_data;
    }
  
    T& data()
    {
      return m_data;
    }

    operator T const &() const
    {
      return m_data;
    }

    operator T&()
    {
      return m_data;
    }
} ;

} // nammespace boost


int main()
{
  Delaunay_triangulation T;

  typedef std::map<Point_2, Vector_3, K::Less_xy_2>           Coord_map;
  typedef CGAL::Data_access<Coord_map>                      Value_access;

  Coord_map value_function;
  double a(0.25), bx(1.3), by(-0.7);

  for (int y=0 ; y<255 ; y++){
    for (int x=0 ; x<255 ; x++){
      K::Point_2_2 p(x,y);
      T.insert(p);
      value_function.insert(std::make_pair(p, Vector_3(x,y,1)));
    }
  }

  //coordinate computation
  K::Point_2_2 p(1.3, 0.34);
  std::vector<std::pair<Point_2, double> > coords;

  double norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
  Vector_3 res =  CGAL::linear_interpolation(coords.begin(), coords.end(), norm,
                                           Value_access(value_function));

  std::cout << "Tested interpolation on " << p << " interpolation: "
            << res << std::endl;

  return EXIT_SUCCESS;
}
