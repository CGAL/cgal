// hds_prog_graph2.C
// -------------------------------------------------
#include <CGAL/Halfedge_data_structure_bases.h>
#include <CGAL/Halfedge_data_structure_using_list.h>
#include <CGAL/Halfedge_data_structure_decorator.h>

using namespace CGAL;

class My_halfedge : public Halfedge_min_base {
    void* prv;
public:
    typedef  Tag_true   Supports_halfedge_prev;
    void*       prev()             { return prv;}
    const void* prev() const       { return prv;}
    void        set_prev( void* h) { prv = h;}
};  

typedef Halfedge_data_structure_using_list <
            Vertex_min_base, My_halfedge, Facet_min_base> HDS;
typedef Halfedge_data_structure_decorator<HDS>  Decorator;

int main() {
    HDS hds;
    Decorator decorator;
    decorator.create_loop( hds);
    CGAL_assertion( decorator.is_valid( hds));
    return 0;
}
