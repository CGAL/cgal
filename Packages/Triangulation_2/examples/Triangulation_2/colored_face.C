#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
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

typedef  CGAL::Cartesian<double> Rp;
typedef  CGAL::Triangulation_euclidean_traits_2<Rp> Gt;
typedef  CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef My_face_base<Gt> Fb;
typedef  CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb > Tds;
typedef  CGAL::Triangulation_2<Gt,Tds> Triangulation;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Face_iterator Face_iterator;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef  CGAL::Point_2<Rp>  Point;

int main() {
    Triangulation t;
    Point p;
   
    while (std::cin >> p){
        t.insert(p);
    }
    Face_iterator fc = t.faces_begin();
    while (fc != t.faces_end()) {
	fc->color = CGAL::BLUE;
	++fc;
    }
    return 0;
}
