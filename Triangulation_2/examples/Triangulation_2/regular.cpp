#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K> Regular_triangulation;

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
