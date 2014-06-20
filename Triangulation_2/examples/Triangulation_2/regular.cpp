/*#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_2<K>  Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;

int main()
{
   std::ifstream in("data/regular.cin");

   Regular_triangulation::Weighted_point wp;
   int count = 0;
   std::vector<Regular_triangulation::Weighted_point> wpoints;
   while(in >> wp){
       count++;
     wpoints.push_back(wp);
   }
   Regular_triangulation rt(wpoints.begin(), wpoints.end());
   rt.is_valid();
   std::cout << "number of inserted points : " << count << std::endl;
   std::cout << "number of vertices :  " ;
   std::cout << rt.number_of_vertices() << std::endl;
   std::cout << "number of hidden vertices :  " ;
   std::cout << rt.number_of_hidden_vertices() << std::endl;
   return 0;
}

*/


#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
//#include <tr1/array>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
int main()
{
  Point_2 points[5] = { Point_2(0,0), Point_2(10,0), Point_2(10,10), Point_2(6,5), Point_2(4,1) };
  Point_2 result[5];
  Point_2 *ptr = CGAL::convex_hull_2( points, points+5, result );
  std::cout <<  ptr - result << " points on the convex hull" << std::endl;
  return 0;
}
