// examples/Halfedge_DS/hds_prog_color.C
// -------------------------------------
#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Halfedge_data_structure_bases.h>
#include <CGAL/Halfedge_data_structure_using_list.h>

using namespace CGAL;

// A facet with a color member variable.
struct My_facet : public Facet_max_base {
    Color color;
};

typedef int  Point;
typedef Halfedge_data_structure_using_list <
    Vertex_max_base<Point>, Halfedge_max_base, My_facet> HDS;

int main() {
    HDS hds;
    My_facet* f = hds.new_facet();
    f->color = RED;
    return 0;
}

    
