// original version see examples/Polyhedron/polyhedron_prog_incr_builder.C
// this time we test compilation with
// an own archetype that uses functionality
// provided by the kernel archetype

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Kernel_archetype.h>

// build an own archetype class ...
// provide 3d point and plane types ...

struct My_archetype {
  typedef CGAL::Point_3_archetype                                   Point_3;
  typedef CGAL::Plane_3_archetype                                   Plane_3;
  typedef CGAL::CGALca::Construct_opposite_plane_3<My_archetype>    Construct_opposite_plane_3;
  
  Construct_opposite_plane_3
  construct_opposite_plane_3_object()
  { return Construct_opposite_plane_3(); }
};


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

typedef My_archetype                       Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

int main() {
    Polyhedron P;
    Build_triangle<HalfedgeDS> triangle;
    P.delegate( triangle);
    CGAL_assertion( P.is_triangle( P.halfedges_begin()));
    return 0;
}
