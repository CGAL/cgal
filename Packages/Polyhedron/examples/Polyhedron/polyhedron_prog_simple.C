// examples/Polyhedron/polyhedron_prog_simple.C
// --------------------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Cartesian<double>            Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    if ( P.is_tetrahedron(h))
        return 0;
    return 1;
}
