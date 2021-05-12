#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/property_map.h>
#include <vector>
#include <cassert>
#include <list>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Partition_traits_2<K, CGAL::Pointer_property_map<K::Point_2>::type > Partition_traits_2;

typedef Partition_traits_2::Point_2                         Point_2;

typedef Partition_traits_2::Polygon_2                       Polygon_2;  // a polygon of indices
typedef std::list<Polygon_2>                                Polygon_list;

/*

      v4     v2
      | \   /|
      |  \ / |
      |  v3  |
      |      |
      v0-----v1

 */

int main( )
{
  std::vector<K::Point_2> points = { K::Point_2(0,0), K::Point_2(2,0), K::Point_2(2,2), K::Point_2(1,1), K::Point_2(0,2) };
  Partition_traits_2 traits(CGAL::make_property_map(points));


  Polygon_2 polygon;
  polygon.push_back(0);
  polygon.push_back(1);
  polygon.push_back(2);
  polygon.push_back(3);
  polygon.push_back(4);

  Polygon_list partition_polys;

  CGAL::y_monotone_partition_2(polygon.vertices_begin(),
                               polygon.vertices_end(),
                               std::back_inserter(partition_polys),
                               traits);

   for (const Polygon_2& poly : partition_polys){
     for(Point_2 p : poly.container()){
        std::cout << "points[" << p << "] =  " << points[p] << ", ";
      }
     std::cout << std::endl;
   }

   assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end(),
                                     traits));

   return 0;
}
