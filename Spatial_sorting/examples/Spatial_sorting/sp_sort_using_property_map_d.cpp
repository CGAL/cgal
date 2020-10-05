#include <CGAL/Cartesian_d.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <vector>

typedef CGAL::Cartesian_d<double>           Kernel;
typedef Kernel::Point_d                     Point_d;
typedef std::pair<Point_d,int>              Point_with_info;
typedef std::vector< Point_with_info >      Data_vector;

//property map and get as friend
// to be allowed to use private member
class Vect_ppmap{
  const Data_vector& points;
public:
  //classical typedefs
  typedef Data_vector::size_type key_type;
  typedef Point_d value_type;
  typedef const value_type& reference;
  typedef boost::readable_property_map_tag category;

  Vect_ppmap(const Data_vector& points_):points(points_){}

  friend reference get(const Vect_ppmap& vmap, key_type i) {
    return vmap.points[i].first;
  }
};

typedef CGAL::Spatial_sort_traits_adapter_d<Kernel,Vect_ppmap>   Search_traits_d;

int main()
{
  double coords[] ={ 1.0, 1.0, 1.0, 1.0,
                     2.0, 2.0, 2.0, 2.0 };

  Data_vector points;
  points.push_back(std::make_pair(Point_d(4,coords  ,coords+4) , 1));
  points.push_back(std::make_pair(Point_d(4,coords+4,coords+8) , 2));

  std::vector<Vect_ppmap::key_type> indices;
  indices.reserve(points.size());

  std::copy(
    boost::counting_iterator<Vect_ppmap::key_type>(0),
    boost::counting_iterator<Vect_ppmap::key_type>(points.size()),
    std::back_inserter(indices) );

  CGAL::spatial_sort(
    indices.begin(),
    indices.end(),
    Search_traits_d(Vect_ppmap(points)) );

  std::vector<Vect_ppmap::key_type>::iterator it=indices.begin();
  for (;it!=indices.end();++it)
    std::cout << points[*it].second << " ";
  std::cout << std::endl;

  std::cout << "done" << std::endl;

  return 0;
}
