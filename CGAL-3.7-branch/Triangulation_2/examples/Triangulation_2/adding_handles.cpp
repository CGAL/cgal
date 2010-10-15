#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <cassert>

/* A vertex class with an additionnal handle */
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class My_vertex_base
  : public  Vb
{
  typedef Vb                              Base;
public:
  typedef typename Vb::Vertex_handle      Vertex_handle;
  typedef typename Vb::Face_handle        Face_handle;
  typedef typename Vb::Point              Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef My_vertex_base<Gt,Vb2>                           Other;
  };

private:
  Vertex_handle  va_;

public:
  My_vertex_base() : Base() {}
  My_vertex_base(const Point & p) : Base(p) {}
  My_vertex_base(const Point & p, Face_handle f) : Base(f,p) {}
  My_vertex_base(Face_handle f) : Base(f) {}

  void set_associated_vertex(Vertex_handle va) { va_ = va;}
  Vertex_handle get_associated_vertex() {return va_ ; }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef My_vertex_base<K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Triangulation_2<K,Tds> Triangulation;

typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Finite_faces_iterator Finite_faces_iterator;
typedef Triangulation::Point   Point;

int main() {
  Triangulation t;
  Vertex_handle v0 = t.insert(Point(0,1));
  Vertex_handle v1 = t.insert(Point(0,0));
  Vertex_handle v2 = t.insert(Point(2,0));
  Vertex_handle v3 = t.insert(Point(2,2));

  // associate vertices as you like
  v0->set_associated_vertex(v1);
  v1->set_associated_vertex(v2);
  v2->set_associated_vertex(v3);
  v3->set_associated_vertex(v0);
  assert( v0->get_associated_vertex() == v1);

  return 0;
}
