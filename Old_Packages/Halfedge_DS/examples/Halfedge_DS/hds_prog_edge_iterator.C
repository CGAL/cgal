// examples/Halfedge_DS/hds_prog_edge_iterator.C
// ---------------------------------------------
#include <CGAL/Halfedge_data_structure_default.h>
#include <CGAL/Halfedge_data_structure_decorator.h>
#include <CGAL/N_step_adaptor.h>

using namespace CGAL;

typedef int                                     Point;
typedef Halfedge_data_structure_default<Point>  HDS;
typedef Halfedge_data_structure_decorator<HDS>  Decorator;
typedef HDS::Halfedge                           Halfedge;
typedef HDS::Halfedge_iterator                  Halfedge_iterator;
typedef HDS::iterator_category                  Iterator_category;
typedef N_step_adaptor< Halfedge_iterator, 2,
                        Halfedge&, Halfedge*,
                        Halfedge,  std::ptrdiff_t,
                        Iterator_category>      Edge_iterator;

int main() {
    HDS hds;
    Decorator decorator;
    decorator.create_loop( hds);
    decorator.create_segment( hds);
    decorator.create_loop( hds);
    CGAL_assertion( decorator.is_valid( hds));
    int n = 0;
    Edge_iterator begin = hds.halfedges_begin();
    for( ; begin != hds.halfedges_end(); ++begin)
        ++n;
    CGAL_assertion( n == 3);  // == 3 edges
    return 0;
}
