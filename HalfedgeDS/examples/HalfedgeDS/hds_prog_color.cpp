#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/IO/Color.h>
#include <cassert>

// A face type with a color member variable.
template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs> {
    CGAL::IO::Color color;
    My_face() {}
    My_face( CGAL::IO::Color c) : color(c) {}
};

// An items type using my face.
struct My_items : public CGAL::HalfedgeDS_items_2 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs> Face;
    };
};

struct My_traits { // arbitrary point type, not used here.
    typedef int  Point_2;
};

typedef CGAL::HalfedgeDS_default<My_traits, My_items> HDS;
typedef HDS::Face                                     Face;
typedef HDS::Face_handle                              Face_handle;

int main() {
    HDS hds;
    Face_handle f = hds.faces_push_back( Face( CGAL::IO::red()));
    f->color = CGAL::IO::blue();
    assert( f->color == CGAL::IO::blue());
    return 0;
}
