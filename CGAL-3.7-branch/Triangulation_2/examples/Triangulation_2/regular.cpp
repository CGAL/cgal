#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_2<K>  Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;

int main()
{
   Regular_triangulation rt;
   std::ifstream in("data/regular.cin");

   Regular_triangulation::Weighted_point wp;
   int count = 0;
   while(in >> wp){
       count++;
     rt.insert(wp);
   }
   rt.is_valid();
   std::cout << "number of inserted points : " << count << std::endl;
   std::cout << "number of vertices :  " ;
   std::cout << rt.number_of_vertices() << std::endl;
   std::cout << "number of hidden vertices :  " ;
   std::cout << rt.number_of_hidden_vertices() << std::endl;
   return 0;
}
