//  file : example/Triangulation_2/triangulation_prog1.C

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

typedef CGAL::Cartesian<double> Rp;
typedef CGAL::Triangulation_euclidean_traits_2<Rp> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt, Tds> Triangulation;
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Gt::Point_2   Point;

int main() {
  std::ifstream in("data/triangulation_prog1.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Triangulation t;
  t.insert(begin, end);
    
  Vertex_circulator vc = t.incident_vertices(t.infinite_vertex()),
    done(vc);
  if (vc != 0) {
    do{
      std::cout << vc->point() << std::endl;
    }while(++vc != done);
  }
  return 0;
}
