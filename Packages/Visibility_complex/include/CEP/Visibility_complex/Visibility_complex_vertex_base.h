#ifndef VISIBILITY_COMPLEX_VERTEX_BASE_H
#define VISIBILITY_COMPLEX_VERTEX_BASE_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef VISIBILITY_COMPLEX_CCW_CW_TRAITS_H
#include <CEP/Visibility_complex/Visibility_complex_ccw_cw_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

template<class _Vc>
class Visibility_complex_vertex_base;

// ------------- Visibility_complex_vertex_base class --------------------------

template<class _Vc>
class Visibility_complex_vertex_base 
    : public _Vc::Gt::Bitangent_2
{
public:
    // -------------------------------------------------------------------------
    typedef typename _Vc::Gt                   Gt;
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    typedef typename Gt::Point_2               Point_2;
    typedef typename Gt::Bitangent_2           Bitangent_2;
    typedef Bitangent_2                        BT;
    typedef typename Gt::Arc_2                 Arc_2;
    typedef typename Gt::Disk                  Disk;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
  //typedef typename Bitangent_2::Disk_handle Disk_handle;
    typedef typename BT::Disk_handle           Disk_handle;
    typedef typename BT::Type                  Type;
    typedef typename BT::Type_util             Type_util;
    // -------------------------------------------------------------------------
    typedef typename _Vc::Vertex               Vertex;
    typedef typename _Vc::Vertex_handle        Vertex_handle;
    typedef typename _Vc::Edge                 Edge;
    typedef typename _Vc::Edge_handle          Edge_handle;
    typedef typename _Vc::Face                 Face;
    typedef typename _Vc::Face_handle          Face_handle;
    // -------------------------------------------------------------------------
public: // private:
    // -------------------------------------------------------------------------
    Edge_handle              _edge[6];
    Vertex_handle            _pi;
    Face_handle              _sup;
    Face_handle              _inf;
    // -------------------------------------------------------------------------
public:    
    // -------------------------------------------------------------------------
    // CONSTRUCTORS
    Visibility_complex_vertex_base() 
	: BT() , 
	  _pi(0) , _sup(0) , _inf(0)
    { 
	for (int i=0;i<6;i++) _edge[i] = 0; 
    }

    Visibility_complex_vertex_base(Type t , Disk_handle start , 
					    Disk_handle finish)
	: BT(t,start,finish) ,
	  _pi(0) , _sup(0) , _inf(0)
    {
	for (int i=0;i<6;i++) _edge[i] = 0; 
    }
    Visibility_complex_vertex_base(Edge_handle start , Edge_handle finish)
    {
	*this = Bitangent_2(Type_util()(start->sign(),finish->sign()),
			    start->object(),finish->object(),
			    *start,*finish);
	for (int i=0;i<6;i++) _edge[i] = 0; 
	_pi = 0; _sup = 0; _inf = 0; 
    }

    Visibility_complex_vertex_base(const Bitangent_2& b) 
	: Bitangent_2(b)
    { 
	for (int i=0;i<6;i++) _edge[i] = 0; 
	_pi = 0; _sup = 0; _inf = 0;
    }
    // -------------------------------------------------------------------------
    // DESTRUCTOR
    ~Visibility_complex_vertex_base() { 
	if (_pi != 0 && 
_pi->_pi == this) 
	  { _pi->_pi = _pi; // test
	  _pi->_pi = 0; }
	if (_sup != 0 && _sup->inf() == this) _sup->set_inf(0);
	if (_inf != 0 && _inf->sup() == this) _inf->set_sup(0);
	for (int i = 0 ; i < 6 ; i++) {
	    if (_edge[i] != 0 && _edge[i]->inf() == this) _edge[i]->set_inf(0);
	    if (_edge[i] != 0 && _edge[i]->sup() == this) _edge[i]->set_sup(0);
	}
    } 
    // ---------- obstacle bitangents ------------------------------------------
    bool is_constraint() const { return (_edge[4] != 0 && _edge[5] != 0); }
    // -------- return the opposite oriented bitangent -------------------------
    Vertex_handle pi();
    void set_pi(Vertex_handle v);
    // -------------------------------------------------------------------------
    // -------- return the four adjacent vertices in ---------------------------
    // ---------- the visiblity graph ------------------------------------------
    Vertex_handle cc(int sens);
    Vertex_handle cc(Disk_handle p , Orientation o = COUNTERCLOCKWISE); 
    // -------------------------------------------------------------------------
    Edge_handle ccw_target_edge() const { return _edge[0]; }
    Edge_handle ccw_source_edge() const { return _edge[1]; }
    Edge_handle ccw_edge(Disk_handle p) const;
    Edge_handle cw_target_edge()  const { return _edge[2]; }
    Edge_handle cw_source_edge()  const { return _edge[3]; }
    Edge_handle cw_edge(Disk_handle p)  const;
    void set_ccw_edge(const Edge_handle& e);
    void set_cw_edge (const Edge_handle& e);
    void set_ccw_edge(const Edge_handle& e , const Disk_handle& p);
    void set_cw_edge (const Edge_handle& e , const Disk_handle& p);
    // -------------------------------------------------------------------------
    Vertex_handle ccR() const 
    { return (ccw_target_edge() == 0) ? 0 : ccw_target_edge()->sup(); }
    Vertex_handle ccL() const 
    { return (ccw_source_edge() == 0) ? 0 : ccw_source_edge()->sup(); }
    Vertex_handle cwR() const 
    { return (cw_source_edge()  == 0) ? 0 :  cw_source_edge()->inf();  }
    Vertex_handle cwL() const 
    { return (cw_target_edge()  == 0) ? 0 : cw_target_edge()->inf();  }
    // -------------------------------------------------------------------------
    Edge_handle target_cusp_edge()  const { return _edge[4]; }
    Edge_handle source_cusp_edge()  const { return _edge[5]; }
    Edge_handle cusp_edge(Disk_handle p) {
	return (p == target_object()) ? _edge[4] : _edge[5];
    }
    void set_target_cusp_edge(Edge_handle e) { _edge[4] = e; }
    void set_source_cusp_edge(Edge_handle e) { _edge[5] = e; }
    void set_cusp_edge(Edge_handle e, Disk_handle p) {
	if (p == target_object()) _edge[4] = e;
	else _edge[5] = e;
    }
    // -------------------------------------------------------------------------
    Face_handle inf()        const { return _inf; }
    void        set_inf(Face_handle f) { _inf = f; }
    Face_handle sup() const { return _sup; }
    void        set_sup(Face_handle f) { _sup = f; }
    // -------------------------------------------------------------------------
    // The two degenerate faces with sink this
    Face_handle target_cusp_face() const 
    {
	if (target_cusp_edge() == 0) return 0;
	return (is_xx_left()) ? target_cusp_edge()->dl(): 
				target_cusp_edge()->dr(); 
    }
    Face_handle source_cusp_face() const
    { 
	if (source_cusp_edge() == 0) return 0;
	return (is_left_xx()) ? source_cusp_edge()->ul() : 
				source_cusp_edge()->ur(); 
    }
    Face_handle cusp_face(Disk_handle p) const
    {
	CGAL_precondition(p == source_object() || p == target_object());
	return (p == source_object()) ? source_cusp_face() : target_cusp_face();
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template< class _Vc >
inline typename Visibility_complex_vertex_base<_Vc>::Edge_handle 
Visibility_complex_vertex_base<_Vc>::ccw_edge(Disk_handle p) const
{
    CGAL_precondition(p == source_object() || p == target_object());
    return (p == source_object()) ? ccw_source_edge() : ccw_target_edge();
}

template< class _Vc >
inline  typename Visibility_complex_vertex_base<_Vc>::Edge_handle 
Visibility_complex_vertex_base<_Vc>::cw_edge(Disk_handle p) const
{
    CGAL_precondition(p == source_object() || p == target_object());
    return (p == source_object()) ? cw_source_edge() : cw_target_edge();
}

template< class _Vc >
inline void
Visibility_complex_vertex_base<_Vc>::set_ccw_edge(const Edge_handle& e , 
						  const Disk_handle& p)
{
    if (target_object() == p) _edge[0] = e;
    else                      _edge[1] = e;
}

template< class _Vc >
inline void
Visibility_complex_vertex_base<_Vc>::set_cw_edge(const Edge_handle& e , 
						 const Disk_handle& p)
{
    if (target_object() == p) _edge[2] = e;
    else                      _edge[3] = e;
}

template< class _Vc >
inline void
Visibility_complex_vertex_base<_Vc>::set_ccw_edge(const Edge_handle& e)
{
    CGAL_precondition(e == 0 || e->object() == source_object() || 
				e->object() == target_object());
    if (e == 0) { _edge[0] = 0; _edge[1] = 0; }
    else if (e->object() == target_object()) _edge[0] = e;
    else _edge[1] = e;
}

template< class _Vc >
inline void
Visibility_complex_vertex_base<_Vc>::set_cw_edge(const Edge_handle& e)
{
    CGAL_precondition(e == 0 || e->object() == source_object() || 
			        e->object() == target_object());
    if (e == 0) { _edge[2] = 0; _edge[3] = 0; }
    else if (e->object() == target_object()) _edge[2] = e;
    else _edge[3] = e;
}

// --------------------------------------------------------------------------------

template< class _Vc >
typename Visibility_complex_vertex_base<_Vc>::Vertex_handle
Visibility_complex_vertex_base<_Vc>::pi()
{
    if (_pi == 0) {
	Type t;
	switch (type()) {
	    case LL: t = RR; break;
	    case LR: t = LR; break;
	    case RL: t = RL; break;
	    case RR: t = LL; break;
	    default: t = LL; break;
	}
	_pi = new Vertex(t ,target_object() , source_object()); 
	_pi->_pi = static_cast<Vertex_handle>(this);
    }
    CGAL_precondition(_pi->_pi == this);
    return _pi;
}

template< class _Vc >
void
Visibility_complex_vertex_base<_Vc>::set_pi(Vertex_handle v)
{
    CGAL_precondition(v != 0);
    if (_pi != 0 && _pi != v) {
	if (_pi->cw_target_edge() != 0) v->set_cw_edge (_pi->cw_target_edge());
	if (_pi->cw_source_edge() != 0) v->set_cw_edge (_pi->cw_source_edge());
	if (_pi->ccw_target_edge() != 0)v->set_ccw_edge(_pi->ccw_target_edge());
	if (_pi->ccw_source_edge() != 0)v->set_ccw_edge(_pi->ccw_source_edge());
	if (_pi->cwL() != 0) _pi->cw_target_edge()->set_sup(v);
	if (_pi->cwR() != 0) _pi->cw_source_edge()->set_sup(v);
	if (_pi->ccR() != 0) _pi->ccw_target_edge()->set_inf(v);
	if (_pi->ccL() != 0) _pi->ccw_source_edge()->set_inf(v);
	if (_pi->inf() != 0) {
	    v->set_inf(_pi->inf());
	    _pi->inf()->set_sup(v);
	}
	if (_pi->sup() != 0) {
	    v->set_sup(_pi->sup());
	    _pi->sup()->set_inf(v);
	}
	delete _pi; 
    }
    _pi = v;
    v->_pi = static_cast<Vertex_handle>(this);
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
