// polyhedron_prog_color.C
// ------------------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

// A face type with a color member variable.
template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs> {
    CGAL::Color color;
};

// An items type using my face.
struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs> Face;
    };
};

typedef CGAL::Cartesian<double>               R;
typedef CGAL::Polyhedron_default_traits_3<R>  Traits;
typedef CGAL::Polyhedron_3<Traits, My_items>  Polyhedron;
typedef Polyhedron::Halfedge_handle           Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    h->facet()->color = CGAL::RED;
    return 0;
}
