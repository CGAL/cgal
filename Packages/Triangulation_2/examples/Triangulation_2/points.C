
#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include "points.h"

#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>


typedef Euclidean_2 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Delaunay_triangulation;

int main(int argc, char* argv[])
{

  Triangulation T;

  PVector V;

  std::ifstream data(argv[1]);
  
  if(data.bad()){
    std::cout << "Problem with file " << argv[1] << std::endl;
  }
  data >> V;

  std::cout << V.size() << std::endl;

  int count = 0;
  std::cout << "Start insertion" << std::endl;
  for(int i = 0; i<V.size();i++){
    
    T.insert(V[i]);

    if(++count == 100){
      std::cout << ".";
      count = 0;
    }
  }
  std::cout << std::endl << "done" << std::endl;
  T.is_valid();

  return 0;
}
