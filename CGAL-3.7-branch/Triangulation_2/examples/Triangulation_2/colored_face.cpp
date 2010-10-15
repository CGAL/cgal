#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<CGAL::Color,K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Triangulation_2<K,Tds> Triangulation;

typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Finite_faces_iterator Finite_faces_iterator;
typedef Triangulation::Point  Point;

int main() {
  Triangulation t;
  t.insert(Point(0,1));
  t.insert(Point(0,0));
  t.insert(Point(2,0));
  t.insert(Point(2,2));

  Finite_faces_iterator fc = t.finite_faces_begin();
  for( ; fc != t.finite_faces_end(); ++fc)  fc->info() = CGAL::BLUE;

  Point p(0.5,0.5);
  Face_handle fh = t.locate(p);
  fh->info() = CGAL::RED;

  return 0;
}
