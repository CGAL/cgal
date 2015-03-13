#include <iostream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Advancing_front_surface_reconstruction<K> Reconstruction;
typedef Reconstruction::Triangulation_3 Triangulation_3;
typedef Reconstruction::TDS_2 TDS_2;
typedef K::Point_3 Point_3;

int main()
{

  std::istream_iterator<Point_3> begin(std::cin);
  std::istream_iterator<Point_3> end;

  Triangulation_3 dt(begin,end);
  
  Reconstruction reconstruction(dt);

  reconstruction.run();
                
  const TDS_2& tds = reconstruction.tds_2();

  for(TDS_2::Face_iterator fit = tds.faces_begin();
      fit != tds.faces_end();
      ++fit){
    if(fit->is_on_surface()){
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
