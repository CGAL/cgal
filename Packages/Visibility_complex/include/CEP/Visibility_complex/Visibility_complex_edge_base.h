#ifndef VISIBILITY_COMPLEX_EDGE_BASE_H
#define VISIBILITY_COMPLEX_EDGE_BASE_H

#include <cmath>
#include <map>

CGAL_BEGIN_NAMESPACE

template< class _Vc >
class Visibility_complex_edge_base;

//------------------------- Visibility_complex_edge_base class ----------------------

template < class _Vc >
class Visibility_complex_edge_base 
    : public _Vc::Gt::Arc_2 
{
public:
    // -------------------------------------------------------------------------
    // Geometry types
    typedef typename _Vc::Gt                    Gt;
    typedef typename Gt::Point_2                Point_2;
    typedef typename Gt::Segment_2              Segment_2;
    typedef typename Gt::Arc_2                  Arc_2;
    typedef Arc_2                               CA;
    typedef typename Gt::Disk                   Disk;
    typedef typename CA::Disk_handle            Disk_handle;
    // -------------------------------------------------------------------------
    typedef typename _Vc::Vertex                Vertex;
    typedef typename _Vc::Vertex_handle         Vertex_handle;
    typedef typename _Vc::Edge                  Edge;
    typedef typename _Vc::Edge_handle           Edge_handle;
    typedef typename _Vc::Face                  Face;
    typedef typename _Vc::Face_handle           Face_handle;
    typedef Visibility_complex_edge_base<_Vc>   Self;
    // -------------------------------------------------------------------------
private:
    bool            _sign;
    Vertex_handle   _sup;
    Vertex_handle   _inf;
public:    Edge_handle     _opposite;
private:    Face_handle     face[3];

public:
    // CONSTRUCTORS ------------------------------------------------------------
    Visibility_complex_edge_base() 
	: CA()   , _sign(true) ,
	  _sup(0) , _inf(0), _opposite(0) 
    { face[0] = 0; face[1] = 0; face[2] = 0; }

    Visibility_complex_edge_base(bool s, Disk_handle p);
    Visibility_complex_edge_base(Vertex_handle v0 , Vertex_handle v1 , 
				 Disk_handle p);
    ~Visibility_complex_edge_base() {
	if (inf() != 0 && this == inf()->source_cusp_edge())
	    inf()->set_source_cusp_edge(0);
	if (inf() != 0 && this == inf()->target_cusp_edge())
	    inf()->set_target_cusp_edge(0);
	if (sup() != 0 && this == sup()->source_cusp_edge())
	    sup()->set_source_cusp_edge(0);
	if (sup() != 0 && this == sup()->target_cusp_edge())
	    sup()->set_target_cusp_edge(0);
	if (sup() != 0 && (object() == sup()->source_object() || 
			   object() == sup()->target_object()) 
		       && this == sup()->cw_edge(object()))
	    sup()->set_cw_edge(0,object());
	if (inf() != 0 && (object() == inf()->source_object() || 
			   object() == inf()->target_object()) 
		       && this == inf()->ccw_edge(object()))
	    inf()->set_ccw_edge(0,object());
	if (face[0] != 0 && this == face[0]->top_edge())
	    face[0]->set_top_edge(0);
	if (face[1] != 0 && this == face[1]->top_edge())
	    face[1]->set_top_edge(0);
	if (face[1] != 0 && this == face[1]->bottom_edge())
	    face[1]->set_bottom_edge(0);
	if (face[2] != 0 && this == face[2]->bottom_edge())
	    face[2]->set_bottom_edge(0);
    }
    // -------------------------------------------------------------------------
    Self& operator=(const Self& e) { Arc_2::operator=(e); return *this; }
    // -------------------------------------------------------------------------
    Face_handle dl() const { return face[0]; }
    Face_handle ul() const { return face[2]; }
    Face_handle dr() const { return (sign()) ? face[0] : face[1]; }
    Face_handle ur() const { return (sign()) ? face[1] : face[2]; }
    void set_adjacent_faces(const Face_handle& f0 , 
			    const Face_handle& f1 , 
			    const Face_handle& f2);
    // -------------------------------------------------------------------------
    bool sign()      const { return _sign;     }
    void set_sign(bool b)  { _sign = b;        }
    void flip_sign()       { _sign = !_sign;   }
    // -------------------------------------------------------------------------
    Vertex_handle  sup()      const { return _sup; }
    void           set_sup(const Vertex_handle& v) { _sup = v; }
    Vertex_handle  inf()      const { return _inf; }
    void           set_inf(const Vertex_handle& v) { _inf = v; }
    // -------------------------------------------------------------------------
    Edge_handle   opposite()  const { return _opposite; }
    void          set_opposite(Edge_handle v);
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template< class _Vc >
Visibility_complex_edge_base<_Vc>::
Visibility_complex_edge_base(bool s, Disk_handle p) : CA(p)
{
    _inf = 0;
    _sup = 0;
    _opposite = 0;
    face[0] = 0; face[1] = 0; face[2] = 0;
    _sign = s;
}

// -----------------------------------------------------------------------------

template< class _Vc >
Visibility_complex_edge_base<_Vc>::
	    Visibility_complex_edge_base(Vertex_handle b0 , Vertex_handle b1 ,
					 Disk_handle p)
    : Arc_2(p,(b0->source_object() == p) ? b0->source() : b0->target(),
		     (b1->source_object() == p) ? b1->source() : b1->target())
{
    _inf = b0;
    _sup = b1;
    _opposite = 0;
    inf()->set_ccw_edge(static_cast<Edge_handle>(this));
    sup()->set_cw_edge (static_cast<Edge_handle>(this));
    face[0] = 0; face[1] = 0; face[2] = 0;

    _sign = (inf()->source_object() == p) ? 
	inf()->is_left_xx() : inf()->is_xx_left();

}

// -----------------------------------------------------------------------------

template< class _Vc >
inline void
Visibility_complex_edge_base<_Vc>::set_opposite(Edge_handle v) 
{ 
    _opposite = v;
    if (v != 0) v->_opposite = static_cast<Edge_handle>(this);
}

// -----------------------------------------------------------------------------

template< class _Vc >
inline void 
Visibility_complex_edge_base<_Vc>::set_adjacent_faces(const Face_handle& f0 , 
						      const Face_handle& f1 , 
						      const Face_handle& f2)
{
    face[0] = f0; face[1] = f1; face[2] = f2;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
