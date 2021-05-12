// computes the normal vector for facets, assuming they are triangles
// (or at least reasonably planar convex polygons), and the normal
// vector for vertices by accumulating the normal vectors of all
// incident facets. All normal vectors are normalized, which requires
// sqrt computations. Therefore we use the Simple_cartesian<double>
// kernel, which is also a natural choice in graphics. Note that the
// normals computed are therefore not exact.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>

// Two functors to compute the normals:  We assume the
// Simple_cartesian<double> Kernel here and use its global functions.

struct Facet_normal {
    template <class Facet>
    void operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        typename Facet::Normal_3 normal = CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - h->next()->vertex()->point());
        f.normal() = normal / std::sqrt( normal * normal);
    }
};

struct Vertex_normal {
    template <class Vertex>
    void operator()( Vertex& v) {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typedef typename Vertex::Halfedge_around_vertex_const_circulator Circ;
        Circ c = v.vertex_begin();
        Circ d = c;
        CGAL_For_all( c, d) {
            if ( ! c->is_border())
                normal = normal + c->facet()->normal();
        }
        v.normal() = normal / std::sqrt( normal * normal);
    }
};

// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).

template <class Refs, class T, class P, class Norm>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
    Norm  norm;
public:
    My_vertex() {} // repeat mandatory constructors
    My_vertex( const P& pt)
      : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
    {
      //initialization of the normal
      norm= CGAL::NULL_VECTOR;
    }
    typedef Norm Normal_3;
    Normal_3&       normal()       { return norm; }
    const Normal_3& normal() const { return norm; }
};

template <class Refs, class T, class Norm>
class My_facet : public CGAL::HalfedgeDS_face_base<Refs, T> {
    Norm  norm;
public:
    // no constructors to repeat, since only default constructor mandatory
    typedef Norm Normal_3;
    Normal_3&       normal()       { return norm; }
    const Normal_3& normal() const { return norm; }
};

struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef My_vertex<Refs, CGAL::Tag_true, Point, Normal> Vertex;
    };
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef typename Traits::Vector_3 Normal;
        typedef My_facet<Refs, CGAL::Tag_true, Normal> Face;
    };
};

// Tie all types together and a small main function using it.

typedef CGAL::Simple_cartesian<double>                 Kernel;
typedef Kernel::Point_3                                Point_3;
typedef CGAL::Polyhedron_3<Kernel, My_items>           Polyhedron;
typedef Polyhedron::Vertex_iterator                    Vertex_iterator;

int main() {
    Point_3 p( 1, 0, 0);
    Point_3 q( 0, 1, 0);
    Point_3 r( 0, 0, 1);
    Point_3 s( 0, 0, 0);
    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    std::for_each( P.facets_begin(),   P.facets_end(),   Facet_normal());
    std::for_each( P.vertices_begin(), P.vertices_end(), Vertex_normal());
    CGAL::set_pretty_mode( std::cout);
    for ( Vertex_iterator i = P.vertices_begin(); i != P.vertices_end(); ++i)
        std::cout << i->normal() << std::endl;
    return 0;
}
