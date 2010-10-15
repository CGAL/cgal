#include <CGAL/HalfedgeDS_min_items.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>

// An items type using a halfedge with previous-pointer.
struct My_items : public CGAL::HalfedgeDS_min_items {
    template <class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef CGAL::HalfedgeDS_halfedge_base< Refs,
            CGAL::Tag_true,  // previous pointer selected.
            CGAL::Tag_false, // pointer to vertex not selected.
            CGAL::Tag_false  // pointer to face not selected.
        > Halfedge;
    };
};

// no traits needed, argument can be arbitrary dummy.
typedef CGAL::HalfedgeDS_default<int, My_items> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>        Decorator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
