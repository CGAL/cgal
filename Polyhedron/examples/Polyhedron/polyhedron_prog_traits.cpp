#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_traits_3<Kernel>  Traits;
typedef CGAL::Polyhedron_3<Traits>         Polyhedron;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    if ( P.is_tetrahedron(h))
        return 0;
    return 1;
}
