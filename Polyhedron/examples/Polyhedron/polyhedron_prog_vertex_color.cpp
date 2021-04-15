#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>

// A vertex type with a color member variable.
template <class Refs, class Point>
struct My_vertex
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>
{
    CGAL::Color color;
    My_vertex() {} // repeat the required constructors
    My_vertex( const Point& p)
        : CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>(p) {}
};

// An items type using my vertex.
struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point_3;
        typedef My_vertex<Refs,Point_3> Vertex;
    };
};

typedef CGAL::Simple_cartesian<double>        Kernel;
typedef CGAL::Polyhedron_3<Kernel, My_items>  Polyhedron;
typedef Polyhedron::Halfedge_handle           Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    h->vertex()->color = CGAL::red();
    return 0;
}
