// polyhedron_prog_incr_builder.C
// -----------------------------------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Cartesian<double>                                R;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R>  HDS;
typedef CGAL::Polyhedron_default_traits_3<R>                   Traits;
typedef CGAL::Polyhedron_3<Traits,HDS>                         Polyhedron;

using CGAL::Modifier_base;

// A modifier creating a triangle using the incremental builder.
template < class HDS>
class Build_triangle : public Modifier_base<HDS> {
public:
    Build_triangle() {}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        typedef typename HDS::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};

int main() {
    Polyhedron P;
    Build_triangle<HDS> triangle;
    P.delegate( triangle);
    CGAL_assertion( P.is_triangle( P.halfedges_begin()));
    return 0;
}
