// file : example/Triangulation_2/points.C

#include <CGAL/basic.h>
#include <fstream>

#include "points.h"
#include <CGAL/Delaunay_triangulation_2.h>

typedef Euclidean_2 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Delaunay_triangulation;

int main()
{
  Triangulation T;
  PVector V;
  std::ifstream data("data/points.cin");
  if(data.bad()){
    std::cout << "Problem with file " << "data/points.cin" << std::endl;
  }
  data >> V;

  std::cout << V.size() << " points read " << std::endl;
  std::cout << "Start insertion" << std::endl;
  for(int i = 0; i<V.size();i++){
    T.insert(V[i]);
  }
  std::cout << std::endl << "done" << std::endl;
  T.is_valid();
  return 0;
}
