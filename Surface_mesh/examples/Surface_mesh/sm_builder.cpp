#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_incremental_builder.h>
#include <CGAL/Surface_mesh.h>

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle {
public:
    Build_triangle() {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Surface_mesh_incremental_builder<HDS> B( hds, true);
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

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

int main() {
    Mesh mesh;
    Build_triangle<Mesh> triangle;
    triangle(mesh);
    //CGAL_assertion( mesh.is_triangle( * mesh.halfedges_begin()));
    return 0;
}
