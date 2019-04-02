#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Pair_partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <cassert>
#include <list>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_2> Surface_mesh;
typedef CGAL::Pair_partition_traits_2<K>            Traits;
typedef Traits::Point_2                                     Point_2;

typedef Traits::Polygon_2                                   Polygon_2;
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

  Surface_mesh sm;

  Traits traits;

  Polygon_2 polygon(traits);

  polygon.push_back(std::make_pair(K::Point_2(0,0), sm.add_vertex(K::Point_2(0,0))));
  polygon.push_back(std::make_pair(K::Point_2(2,0), sm.add_vertex(K::Point_2(2,0))));
  polygon.push_back(std::make_pair(K::Point_2(2,2), sm.add_vertex(K::Point_2(2,2))));
  polygon.push_back(std::make_pair(K::Point_2(1,1), sm.add_vertex(K::Point_2(1,1))));
  polygon.push_back(std::make_pair(K::Point_2(0,2), sm.add_vertex(K::Point_2(0,2))));
  
  Polygon_list partition_polys;



  //CGAL::y_monotone_partition_2
  CGAL::optimal_convex_partition_2
     (polygon.vertices_begin(),
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys),
                                traits);

   std::list<Polygon_2>::const_iterator   poly_it;
   for (poly_it = partition_polys.begin(); poly_it != partition_polys.end();
        poly_it++)
   {
     for(Point_2 p : poly_it->container()){
        std::cout << "[" << p.first << "| " << p.second << "] ";
      }
     std::cout << std::endl;
   }
#if 1
   assert(CGAL::convex_partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end(),
                                     traits));
#endif
   return 0;
}
