// examples/Halfedge_DS/hds_prog_halfedge_iterator.C
// -------------------------------------------------
#include <CGAL/Halfedge_data_structure_default.h>
#include <CGAL/Halfedge_data_structure_decorator.h>

typedef int                                           Point;
typedef CGAL::Halfedge_data_structure_default<Point>  HDS;
typedef CGAL::Halfedge_data_structure_decorator<HDS>  Decorator;
typedef HDS::Halfedge_iterator                        Halfedge_iterator;

int main() {
    HDS hds;
    Decorator decorator;
    decorator.create_loop( hds);
    decorator.create_segment( hds);
    decorator.create_loop( hds);
    CGAL_assertion( decorator.is_valid( hds));
    int n = 0;
    Halfedge_iterator begin = hds.halfedges_begin();
    for( ; begin != hds.halfedges_end(); ++begin)
        ++n;
    CGAL_assertion( n == 6);  // == 3 edges
    return 0;
}
