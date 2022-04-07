#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

template < class GT, class Vb = CGAL::Triangulation_vertex_base_3<GT> >
class My_vertex_base
  : public Vb
{
public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Cell_handle    Cell_handle;
  typedef typename Vb::Point          Point;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<GT, Vb2>                        Other;
  };

  My_vertex_base() {}
  My_vertex_base(const Point& p) : Vb(p) {}
  My_vertex_base(const Point& p, Cell_handle c) : Vb(p, c) {}

  Vertex_handle   vh;
  Cell_handle     ch;
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_cell_base_3<K>                  Cb;
typedef CGAL::Triangulation_data_structure_3<My_vertex_base<K>, Cb>  Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                       Delaunay;

typedef Delaunay::Vertex_handle                                      Vertex_handle;
typedef Delaunay::Point                                              Point;

int main()
{
  Delaunay T;

  Vertex_handle v0 = T.insert(Point(0,0,0));
  Vertex_handle v1 = T.insert(Point(1,0,0));
  Vertex_handle v2 = T.insert(Point(0,1,0));
  Vertex_handle v3 = T.insert(Point(0,0,1));
  Vertex_handle v4 = T.insert(Point(2,2,2));
  Vertex_handle v5 = T.insert(Point(-1,0,1));

  // Now we can link the vertices as we like.
  v0->vh = v1;
  v1->vh = v2;
  v2->vh = v3;
  v3->vh = v4;
  v4->vh = v5;
  v5->vh = v0;

  return 0;
}
