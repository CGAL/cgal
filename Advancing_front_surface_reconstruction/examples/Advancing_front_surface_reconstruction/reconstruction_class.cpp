#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Advancing_front_surface_reconstruction<K> Reconstruction;
typedef Reconstruction::Triangulation_3 Triangulation_3;
typedef Reconstruction::Triangulation_data_structure_2 TDS_2;
typedef K::Point_3 Point_3;

int main(int argc, char* argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/half.xyz");
  std::istream_iterator<Point_3> begin(in);
  std::istream_iterator<Point_3> end;

  Triangulation_3 dt(begin,end);
  
  Reconstruction reconstruction(dt);

  reconstruction.run();
                
  const TDS_2& tds = reconstruction.triangulation_data_structure_2();

  for(TDS_2::Face_iterator fit = tds.faces_begin();
      fit != tds.faces_end();
      ++fit){
    if(reconstruction.has_on_surface(fit)){
      Triangulation_3::Facet f = fit->facet();
      Triangulation_3::Cell_handle ch = f.first;
      int ci = f.second;
      for(int i = 0; i < 4; i++){
        if(ci != i){
          std:: cout << ch->vertex(i)->point() << "   ";
        }
      }
      std::cout << std::endl;
    }
  }
  
  return 0;
}
