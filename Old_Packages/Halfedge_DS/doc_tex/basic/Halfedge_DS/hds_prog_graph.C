// hds_prog_graph.C
// -------------------------------------------------
#include <CGAL/Halfedge_data_structure_bases.h>
#include <CGAL/Halfedge_data_structure_using_list.h>
#include <CGAL/Halfedge_data_structure_decorator.h>

typedef CGAL::Halfedge_data_structure_using_list <
    CGAL::Vertex_min_base, CGAL::Halfedge_min_base, CGAL::Facet_min_base> HDS;
typedef CGAL::Halfedge_data_structure_decorator<HDS>  Decorator;

int main() {
    HDS hds;
    Decorator decorator;
    decorator.create_loop( hds);
    CGAL_assertion( decorator.is_valid( hds));
    return 0;
}
