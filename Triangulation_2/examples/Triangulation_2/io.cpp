#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <fstream>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

struct Height {
  int height;
};

std::ostream& operator<<(std::ostream& os, const Height& h)
{
  os << h.height;
  return os;
}

std::istream& operator>>(std::istream& is, Height& h)
{
  is >> h.height;
  return is;
}

typedef CGAL::Triangulation_vertex_base_with_info_2<Height,K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<CGAL::Color,K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Triangulation_2<K,Tds> Triangulation;

typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Finite_faces_iterator Finite_faces_iterator;
typedef Triangulation::Point  Point;

int main() {
  Triangulation t;
  t.insert(Point(0,1))->info().height = 78;
  t.insert(Point(0,0))->info().height = 12;
  t.insert(Point(2,0))->info().height = 0;
  t.insert(Point(2,2))->info().height = 7812;

  Finite_faces_iterator fc = t.finite_faces_begin();
  for( ; fc != t.finite_faces_end(); ++fc){
    fc->info() = CGAL::BLUE;
  }
  Point p(0.5,0.5);
  Face_handle fh = t.locate(p);
  fh->info() = CGAL::RED;

  {
    std::ofstream os("triangulation.txt");
    os << t << std::endl;
  }

  t.clear();
  {
    std::ifstream is("triangulation.txt");
    is >> t;
  }

  std::cout << t << std::endl; 
  

  
  return 0;
}
