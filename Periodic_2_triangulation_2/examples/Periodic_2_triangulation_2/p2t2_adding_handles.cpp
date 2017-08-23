#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

template < class GT, class Vb >
class My_vertex_base
  : public Vb
{
public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Face_handle    Face_handle;
  typedef typename Vb::Point          Point;

  template < typename Tds2 >
  struct Rebind_TDS
  {
    typedef typename Vb::template Rebind_TDS<Tds2>::Other     Vb2;
    typedef My_vertex_base<GT, Vb2>                           Other;
  };

  My_vertex_base() {}

  My_vertex_base(const Point& p)
    : Vb(p) {}

  My_vertex_base(const Point& p, Face_handle c)
    : Vb(p, c) {}

  Vertex_handle   vh;
  Face_handle     fh;
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT>    VbDS;
typedef My_vertex_base<GT, VbDS>                            Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2<GT>      Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds>  PDT;

typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Point                                          Point;

int main()
{
  PDT T;

  Vertex_handle v0 = T.insert(Point(0, 0));
  Vertex_handle v1 = T.insert(Point(.1, 0));
  Vertex_handle v2 = T.insert(Point(0, .1));
  Vertex_handle v3 = T.insert(Point(0, 0.2));
  Vertex_handle v4 = T.insert(Point(.2, .2));
  Vertex_handle v5 = T.insert(Point(.9, 0));

  // Now we can link the vertices as we like.
  v0->vh = v1;
  v1->vh = v2;
  v2->vh = v3;
  v3->vh = v4;
  v4->vh = v5;
  v5->vh = v0;

  return 0;
}
