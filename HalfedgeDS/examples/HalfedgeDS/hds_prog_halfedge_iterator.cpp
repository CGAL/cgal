#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <cassert>

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_default<Traits> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;
typedef HDS::Halfedge_iterator           Iterator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    decorator.create_segment();
    assert( decorator.is_valid());
    int n = 0;
    for ( Iterator i = hds.halfedges_begin(); i != hds.halfedges_end(); ++i )
        ++n;
    assert( n == 4);  // == 2 edges
    return 0;
}
