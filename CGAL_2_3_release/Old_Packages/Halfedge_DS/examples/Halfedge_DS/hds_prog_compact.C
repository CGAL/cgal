// examples/Halfedge_DS/hds_prog_compact.C
// ---------------------------------------
#include <CGAL/Halfedge_data_structure_bases.h>
#include <CGAL/Halfedge_data_structure_using_vector.h>
#include <CGAL/Halfedge_data_structure_decorator.h>

using namespace CGAL;

// Define a new halfedge class.
class My_halfedge {
protected:
    size_t  nxt;
    void*   v;
    void*   f;
public:
    typedef Tag_false Supports_halfedge_prev;
    typedef Tag_true  Supports_halfedge_vertex;
    typedef Tag_true  Supports_halfedge_facet;

    My_halfedge() : nxt(0), f(NULL) {}

    void*       opposite()       {
        const size_t SIZE = sizeof( My_halfedge);
        if ( nxt & 1)
            return (char*)this + SIZE;
        return (char*)this - SIZE;
    }
    const void* opposite() const {
        const size_t SIZE = sizeof( My_halfedge);
        if ( nxt & 1)
            return (const char*)this + SIZE;
        return (const char*)this - SIZE;
    }
    void*       next()           { return (void*)(nxt & (~ size_t(1)));}
    const void* next() const     { return (void*)(nxt & (~ size_t(1)));}

    void*       vertex()         { return v;}
    const void* vertex() const   { return v;}

    void*       facet()          { return f;}
    const void* facet() const    { return f;}

    bool is_border() const       { return f == NULL;}

    void  set_opposite( void* g) {
        char* h = (char*)g;
        CGAL_assertion( size_t( abs( h - (char*)this)) == sizeof(My_halfedge));
        if ( h > (char*)this)
            nxt |= 1;
        else
            nxt &= (~ size_t(1));
    }
    void  set_next( void* h)     {
        CGAL_assertion( ((size_t)h & 1) == 0);
        nxt = ((size_t)(h)) | (nxt & 1);
    }
    void  set_vertex( void* _v)  { v = _v;}
    void  set_facet( void* _f)   { f = _f;}
};

typedef int  Point;
typedef Halfedge_data_structure_using_vector <
            Vertex_max_base<Point>, My_halfedge, Facet_max_base> HDS;
typedef Halfedge_data_structure_decorator<HDS>  Decorator;

int main() {
    HDS hds(1,2,2);
    Decorator decorator;
    decorator.create_loop( hds);
    CGAL_assertion( decorator.is_valid(hds));
    return 0;
}
