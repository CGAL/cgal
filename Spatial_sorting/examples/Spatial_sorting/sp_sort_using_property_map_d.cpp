#include <CGAL/Cartesian_d.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <boost/iterator/counting_iterator.hpp>
#include <vector>

typedef CGAL::Cartesian_d<double>           Kernel;
typedef Kernel::Point_d                     Point_d;
typedef std::pair<Point_d,int>              Point_with_info;
typedef std::vector< Point_with_info >      Data_vector;

//property map and get as friend
// to be allowed to use private member
class Vect_ppmap{
  const Data_vector& points
public:
  //classical typedefs
  typedef Data_vector::size_t key_type;
  typedef Point_d value_type;
  typedef const Point_d& reference;
  typedef boost::readable_property_map_tag category;

  Vect_ppmap(const Data_vector& points_):points(points_){}

  friend reference get(const Vect_ppmap& vmap,key_type i) const{
    return vmap.points[i];
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

  std::vector<std::ptrdiff_t> indices;
  indices.reserve(points.size());  
  
  std::copy(boost::counting_iterator<std::ptrdiff_t>(0),
            boost::counting_iterator<std::ptrdiff_t>(points.size()),
            std::back_inserter(indices));
  
  CGAL::spatial_sort( indices.begin(),indices.end(),Search_traits_d(Vect_ppmap(points)) );
  
  for (Data_vector::iterator it=points.begin();it!=points.end();++it)
    std::cout << it->second << " ";
  std::cout << "\n";

  std::cout << "done" << std::endl;
  
  return 0;
}
