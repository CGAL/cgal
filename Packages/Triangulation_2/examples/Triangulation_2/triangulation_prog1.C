#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

using namespace CGAL;

typedef Cartesian<double> Rp;
typedef Triangulation_euclidean_traits_2<Rp> Gt;
typedef Triangulation_vertex_base_2<Gt> Vb;
typedef Triangulation_face_base_2<Gt>  Fb;
typedef Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef Triangulation_2<Gt, Tds> Triangulation;
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Gt::Point_2   Point;

int main() {
  Triangulation t;
    
  Point p;
  while (std::cin >> p){
    t.insert(p);
  }
  
  Vertex_circulator vc = t.incident_vertices(t.infinite_vertex()),
    done(vc);
  if (vc != 0) {
    do{
      std::cout << vc->point() << std::endl;
    }while(++vc != done);
  }
  return 0;
}
