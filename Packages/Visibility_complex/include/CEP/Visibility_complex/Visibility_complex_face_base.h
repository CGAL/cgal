#ifndef VISIBILITY_COMPLEX_FACE_BASE_H
#define VISIBILITY_COMPLEX_FACE_BASE_H

#include <cmath>
//#include <limits.h>
#include <queue>

#ifndef VISIBILITY_COMPLEX_CCW_CW_TRAITS_H
#include <CEP/Visibility_complex/Visibility_complex_ccw_cw_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template< class _Vc >
class Visibility_complex_face_base {
private:
    typedef typename _Vc::Gt                   Gt;
    typedef typename Gt::Disk                  Disk;
    typedef typename _Vc::Disk_handle          Disk_handle;
    typedef typename _Vc::Vertex               Vertex;
    typedef typename _Vc::Vertex_handle        Vertex_handle;
    typedef typename _Vc::Edge                 Edge;
    typedef typename _Vc::Edge_handle          Edge_handle;
    typedef typename _Vc::Face                 Face;
    typedef typename _Vc::Face_handle          Face_handle;

private:
    Vertex_handle   _sup;
    Vertex_handle   _inf;

    Edge_handle _back_view;
    Edge_handle _front_view;
    Disk_handle _back_object;
    Disk_handle _front_object;

    Edge_handle     _bottom_edge;
    Edge_handle     _top_edge;

public:
    // CONSTRUCTORS  -----------------------------------------------------------
    Visibility_complex_face_base() 
	: _sup(0)         , _inf(0),
	  _back_view(0)   , _front_view(0) ,
	  _back_object(0) , _front_object(0) ,
	  _bottom_edge(0) , _top_edge(0)
    { }

    Visibility_complex_face_base(Edge_handle s     , Edge_handle t,
				 Disk_handle ba , Disk_handle fr)
	: _sup(0)         , _inf(0)           ,
	  _back_view(0)   , _front_view(0) ,
	  _back_object(ba), _front_object(fr) ,
	  _bottom_edge(s) , _top_edge(t)
    { }
    // -------------------------------------------------------------------------
    // DESTRUCTOR
    ~Visibility_complex_face_base() { 
	if (sup() != 0 && this == sup()->inf()) sup()->set_inf(0);	
	if (inf() != 0 && this == inf()->sup()) inf()->set_sup(0);	
	Edge_handle e = bottom_edge();
	if (e != 0) {
	    if (!e->sign() && e->ul() == this) 
		e->set_adjacent_faces(e->dl(),e->dr(),0);
	    else if (e->sign() && e->ul() == this)
		e->set_adjacent_faces(e->dl(),e->ur(),0);
	    else if (e->sign() && e->ur() == this)
		e->set_adjacent_faces(e->dl(),0,e->ul());
	}
	Edge_handle f = top_edge();
	if (f != 0) {
	    if (f->sign() && f->dl() == this)
		f->set_adjacent_faces(0,f->ur(),f->ul());
	    else if (!f->sign() && f->dl() == this)
		f->set_adjacent_faces(0,f->dr(),f->ul());
	    else if (!f->sign() && f->dr() == this)
		f->set_adjacent_faces(f->dl(),0,f->ul());
	}
    }
    // -------------------------------------------------------------------------
    Vertex_handle sup() const { return _sup; }
    Vertex_handle inf() const { return _inf; }
    void set_sup(const Vertex_handle& v);
    void set_inf(const Vertex_handle& v);
    // -------------------------------------------------------------------------
    Edge_handle front_view()  const    { return _front_view; }
    Edge_handle back_view()   const    { return _back_view;  }
    void set_back_view (Edge_handle e) { _back_view  = e;    }
    void set_front_view(Edge_handle e) { _front_view = e;    }
    // -------------------------------------------------------------------------
    Disk_handle front_object()  const    { return _front_object; }
    Disk_handle back_object()   const    { return _back_object;  }
    void set_back_object (Disk_handle p) { _back_object  = p;    }
    void set_front_object(Disk_handle p) { _front_object = p;    }
    // -------------------------------------------------------------------------
    void set_top_edge(const Edge_handle& e)    { _top_edge    = e;     }
    void set_bottom_edge(const Edge_handle& e) { _bottom_edge = e;     }
    Edge_handle top_edge()      const          { return _top_edge;     }
    Edge_handle bottom_edge()   const          { return _bottom_edge;  }
    // -------------------------------------------------------------------------
    //template < class _Tr > bool candidate(_Tr) const; 
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc >
void 
Visibility_complex_face_base<_Vc>::set_sup(const Vertex_handle& v)
{
    _sup = v;
    if (v != 0 && !v->is_constraint()) 
	v->set_inf(static_cast<Face_handle>(this));
}

template < class _Vc >
void 
Visibility_complex_face_base<_Vc>::set_inf(const Vertex_handle& v)
{
    _inf = v;
    if (v != 0 && !v->is_constraint()) 
	v->set_sup(static_cast<Face_handle>(this));
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
