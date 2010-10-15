#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/N_step_adaptor.h>

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_default<Traits>            HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>             Decorator;
typedef HDS::Halfedge_iterator                      Halfedge_iterator;
typedef CGAL::N_step_adaptor< Halfedge_iterator, 2> Iterator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    decorator.create_segment();
    CGAL_assertion( decorator.is_valid());
    int n = 0;
    for ( Iterator e = hds.halfedges_begin(); e != hds.halfedges_end(); ++e)
        ++n;
    CGAL_assertion( n == 2);  // == 2 edges
    return 0;
}
