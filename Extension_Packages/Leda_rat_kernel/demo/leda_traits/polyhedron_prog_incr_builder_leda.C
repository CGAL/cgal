// see examples/Polyhedron/polyhedron_prog_incr_builder.C

// provide 3d kernel traits ...
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CGAL/Cartesian.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

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
        B.add_vertex( Point( 0, 0, 0, 1));
        B.add_vertex( Point( 1, 0, 0, 1));
        B.add_vertex( Point( 0, 1, 0, 1));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};

//typedef CGAL::Cartesian<double>            Kernel;
typedef CGAL::leda_rat_kernel_traits       Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

int main() {
    Polyhedron P;
    Build_triangle<HalfedgeDS> triangle;
    P.delegate( triangle);
    CGAL_assertion( P.is_triangle( P.halfedges_begin()));
    return 0;
}
