// file example/Triangulation_2/regular.C
#include <CGAL/Cartesian.h>
#include <fstream>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

typedef CGAL::Cartesian<double> Rp;
typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Rp,W>  Gt;
typedef CGAL::Regular_triangulation_2<Gt> Regular_triangulation;

int main()
{
   Regular_triangulation rt;
   std::ifstream in("data/regular.cin");

   Gt::Weighted_point wp;
   while(in >> wp){
     std::cout << wp << std::endl;
     rt.insert(wp);
   }
   rt.is_valid();
   return 0;	
}
