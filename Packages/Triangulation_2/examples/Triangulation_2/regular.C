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
   while(in >> wp){
     std::cout << wp << std::endl;
     rt.insert(wp);
     rt.is_valid();
   }
   rt.is_valid();
   return 0;	
}
