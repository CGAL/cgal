// file example/Triangulation_2/regular.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <fstream>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<K,W>  Gt;
typedef CGAL::Regular_triangulation_2<Gt> Regular_triangulation;

int main()
{
   Regular_triangulation rt;
   std::ifstream in("data/regular.cin");

   Gt::Weighted_point wp;
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
