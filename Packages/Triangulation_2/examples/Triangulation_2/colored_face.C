// file : examples/Triangulation_2/colored_face.C
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>

/* A facet with a color member variable. */
template < class Fb = CGAL::Triangulation_ds_face_base_2<> >
class My_face_base 
  : public  Fb
{
  typedef Fb                                           Base;
  typedef typename Fb::Triangulation_data_structure    TDS;
public:
  typedef TDS                              Triangulation_data_structure;
  typedef typename TDS::Vertex_handle      Vertex_handle;
  typedef typename TDS::Face_handle        Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other    Fb2;
    typedef My_face_base<Fb2>                                  Other;
  };

  CGAL::Color color;

  My_face_base() : Base() {}

  My_face_base(Vertex_handle v0, 
	       Vertex_handle v1, 
	       Vertex_handle v2) : Base(v0,v1,v2) {}

  My_face_base(Vertex_handle v0, 
	       Vertex_handle v1, 
	       Vertex_handle v2,
	       Face_handle n0, 
	       Face_handle n1, 
	       Face_handle n2) : Base(v0,v1,v2,n0,n1,n2) {}
};

typedef CGAL::Cartesian<double> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef My_face_base<> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb > Tds;
typedef CGAL::Triangulation_2<Gt,Tds> Triangulation;

typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Finite_faces_iterator Finite_faces_iterator;
typedef Gt::Point_2  Point;

int main() {
  Triangulation t;
  t.insert(Point(0,1));
  t.insert(Point(0,0));
  t.insert(Point(2,0));
  t.insert(Point(2,2));
 
  Finite_faces_iterator fc = t.finite_faces_begin();
  for( ; fc != t.finite_faces_end(); ++fc)  fc->color = CGAL::BLUE;

  Point p(0.5,0.5);
  Face_handle fh = t.locate(p);
  fh->color = CGAL::RED;

  return 0;
}
