#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>


typedef CGAL::Cartesian<double> Rp;
typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Rp,W>  Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Regular_triangulation_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef CGAL::Regular_triangulation_2<Gt, Tds> Regular_triangulation;

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
