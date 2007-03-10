// original version see examples/Polyhedron/polyhedron_prog_incr_builder.cpp
// this time we test compilation with the
// kernel concept archetype
// we use interface restrictions in this example

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

// limit interface of the kernel archetype ...
#define CGAL_CA_LIMITED_INTERFACE
#define CGAL_CA_POINT_3
#define CGAL_CA_PLANE_3

#include <CGAL/Kernel_archetype.h>

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
    Build_triangle() {}
    void operator()( HDS& hds) {
        // Postcondition: `hds' is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point());
        B.add_vertex( Point());
        B.add_vertex( Point());
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};

typedef CGAL::Kernel_archetype             Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

int main() {
    Polyhedron P;
    Build_triangle<HalfedgeDS> triangle;
    P.delegate( triangle);
    CGAL_assertion( P.is_triangle( P.halfedges_begin()));
    return 0;
}
