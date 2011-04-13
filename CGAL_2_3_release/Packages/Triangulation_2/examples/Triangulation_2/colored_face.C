// file : examples/Triangulation_2/colored_face.C
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>

/* A facet with a color member variable. */
template < class Gt >
class My_face_base : public CGAL::Triangulation_face_base_2<Gt>
{
public:
  CGAL::Color color;
  My_face_base() :
    CGAL::Triangulation_face_base_2<Gt>() {}
  My_face_base(void* v0, void* v1, void* v2) : 
    CGAL::Triangulation_face_base_2<Gt>(v0,v1,v2) {}
  My_face_base(void* v0, void* v1, void* v2, void* n0, void* n1, void* n2) : 
    CGAL::Triangulation_face_base_2<Gt>(v0,v1,v2,n0,n1,n2) {}
};

typedef CGAL::Cartesian<double> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef My_face_base<Gt> Fb;
typedef CGAL::Triangulation_data_structure_using_list_2<Vb,Fb > Tds;
typedef CGAL::Triangulation_2<Gt,Tds> Triangulation;

typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Face_iterator Face_iterator;
typedef Gt::Point_2  Point;

int main() {
  Triangulation t;
  t.insert(Point(0,1));
  t.insert(Point(0,0));
  t.insert(Point(2,0));
  t.insert(Point(2,2));
 
  Face_iterator fc = t.faces_begin();
  for( ; fc != t.faces_end(); ++fc)  fc->color = CGAL::BLUE;

  Point p(0.5,0.5);
  Face_handle fh = t.locate(p);
  fh->color = CGAL::RED;

  return 0;
}
