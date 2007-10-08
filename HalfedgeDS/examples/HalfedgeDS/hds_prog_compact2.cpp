#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <cstddef>

// Define a new halfedge class. We assume that the Halfedge_handle can
// be created from a pointer (e.g. the HalfedgeDS is based here on the
// In_place_list internally) and that halfedges are allocated in pairs.
// We encode the opposite pointer in a single bit, which is stored
// in the lower bit of the next-pointer. We use the static member function
// HDS::halfedge_handle to translate pointer to handles.
template <class Refs>
class My_halfedge {
public:
    typedef Refs                                 HDS;
    typedef My_halfedge<Refs>                    Base_base;
    typedef My_halfedge<Refs>                    Base;
    typedef My_halfedge<Refs>                    Self;
    typedef CGAL::Tag_false                      Supports_halfedge_prev;
    typedef CGAL::Tag_true                       Supports_halfedge_vertex;
    typedef CGAL::Tag_true                       Supports_halfedge_face;
    typedef typename Refs::Vertex_handle         Vertex_handle;
    typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
    typedef typename Refs::Halfedge              Halfedge;
    typedef typename Refs::Halfedge_handle       Halfedge_handle;
    typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Refs::Face_handle           Face_handle;
    typedef typename Refs::Face_const_handle     Face_const_handle;
private:
    std::ptrdiff_t  nxt;
public:
    My_halfedge() : nxt(0), f( Face_handle()) {}

    Halfedge_handle opposite() {
        // Halfedge could be different from My_halfedge (e.g. pointer for
        // linked list). Get proper handle from 'this' pointer first, do
        // pointer arithmetic, then convert pointer back to handle again.
        Halfedge_handle h = HDS::halfedge_handle(this); // proper handle
        if ( nxt & 1)
            return HDS::halfedge_handle( &* h + 1);
        return HDS::halfedge_handle( &* h - 1);
    }
    Halfedge_const_handle opposite() const { // same as above
        Halfedge_const_handle h = HDS::halfedge_handle(this); // proper handle
        if ( nxt & 1)
            return HDS::halfedge_handle( &* h + 1);
        return HDS::halfedge_handle( &* h - 1);
    }
    Halfedge_handle next() {
        return HDS::halfedge_handle((Halfedge*)(nxt & (~ std::ptrdiff_t(1))));
    }
    Halfedge_const_handle next() const {
        return Halfedge_const_handle((const Halfedge*)(nxt &
                                                 (~ std::ptrdiff_t(1))));
    }
    void  set_opposite( Halfedge_handle h) {
        CGAL_precondition(( &* h - 1 == &* HDS::halfedge_handle(this)) ||
                          ( &* h + 1 == &* HDS::halfedge_handle(this)));
        if ( &* h - 1 == &* HDS::halfedge_handle(this))
            nxt |= 1;
        else
            nxt &= (~ std::ptrdiff_t(1));
    }
    void  set_next( Halfedge_handle h) {
        CGAL_precondition( ((std::ptrdiff_t)(&*h) & 1) == 0);
        nxt = ((std::ptrdiff_t)(&*h)) | (nxt & 1);
    }
private:    // Support for the Vertex_handle.
    Vertex_handle    v;
public:
    // the incident vertex.
    Vertex_handle         vertex()                     { return v; }
    Vertex_const_handle   vertex() const               { return v; }
    void                  set_vertex( Vertex_handle w) { v = w; }

private:
    Face_handle      f;
public:
    Face_handle           face()                       { return f; }
    Face_const_handle     face() const                 { return f; }
    void                  set_face( Face_handle g)     { f = g; }
    bool                  is_border() const { return f == Face_handle(); }
};

// Replace halfedge in the default items type.
struct My_items : public CGAL::HalfedgeDS_items_2 {
    template <class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef My_halfedge<Refs> Halfedge;
    };
};

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_list<Traits, My_items> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
