#ifndef VISIBILITY_COMPLEX_CCW_CW_TRAITS_H
#define VISIBILITY_COMPLEX_CCW_CW_TRAITS_H

#include <functional>
#include <CEP/Visibility_complex/Visibility_complex_function_objects.h>

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------- 

template < class _Vc >
struct Visibility_complex_duality_traits_base
{
    typedef typename _Vc::Gt             Gt;
    typedef typename _Vc::Disk_handle Disk_handle;
    typedef typename _Vc::Vertex         Vertex;
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    typedef typename _Vc::Edge           Edge;
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Face_handle    Face_handle;
    // ------------------------------------------------------------------------- 
    struct Sup {
	typedef Vertex_handle result_type;
	template <class _Arg >
	Vertex_handle operator()(const _Arg& a) const { return a->sup(); }
	Face_handle   operator()(const Vertex_handle& v) const 
	{ return v->sup(); }
    };
    // ------------------------------------------------------------------------- 
    struct Set_sup {
	void operator()(const Vertex_handle& v, Face_handle f) const 
	{ v->set_sup(f); }
	template <class _Arg >
	void operator()(_Arg a,const Vertex_handle& v) const { a->set_sup(v); }
    };
    // ------------------------------------------------------------------------- 
    struct Inf {
	typedef Vertex_handle result_type;
	template <class _Arg >
	Vertex_handle operator()(const _Arg& a)   const { return a->inf(); }
	Face_handle   operator()(const Vertex_handle& v) const 
	{ return v->inf(); }
    };
    // ------------------------------------------------------------------------- 
    struct Set_inf {
	void operator()(Vertex_handle v,Face_handle f) const { v->set_inf(f); }
	template <class _Arg >
	void operator()(_Arg a,Vertex_handle v) const { a->set_inf(v); }
    };
    // ------------------------------------------------------------------------- 
    struct Ccw_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->ccw_edge(p); }
    };
    // -------------------------------------------------------------------------
    struct Cw_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->cw_edge(p); }
    };
    // -------------------------------------------------------------------------
    struct Cusp_edge {
	Edge_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->cusp_edge(p); }
    };
    // -------------------------------------------------------------------------
    struct Set_cusp_edge {
	void operator()(const Vertex_handle& v, 
			Edge_handle e ,const Disk_handle& p) const
	{ v->set_cusp_edge(e,p); }
    };
    // -------------------------------------------------------------------------
    struct Cusp_face {
	Face_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->cusp_face(p); }
    };
    // -------------------------------------------------------------------------
    struct Set_ccw_edge {
	void operator()(const Vertex_handle& v, Edge_handle e) const
	{ v->set_ccw_edge(e); }
    };
    // -------------------------------------------------------------------------
    struct Set_cw_edge {
	void operator()(const Vertex_handle& v, Edge_handle e) const
	{ v->set_cw_edge(e); }
    };
    // -------------------------------------------------------------------------
    struct Source_object 
	: public std::unary_function<Vertex_handle,Disk_handle> {
	Disk_handle operator()(const Vertex_handle& v) const 
	{ return v->source_object(); }
    };
    // -------------------------------------------------------------------------
    struct Target_object
	: public std::unary_function<Vertex_handle,Disk_handle> {
	Disk_handle operator()(const Vertex_handle& v) const 
	{ return v->target_object(); }
    };
    // -------------------------------------------------------------------------
    struct Top_edge
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& v) const 
	{ return v->top_edge(); }
    };
    // -------------------------------------------------------------------------
    struct Bottom_edge
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& v) const 
	{ return v->bottom_edge(); }
    };
    // -------------------------------------------------------------------------
    struct Dl : public std::unary_function<Edge_handle,Face_handle> { 
	Face_handle operator()(const Edge_handle& e) const { return e->dl(); }
    };
    // ------------------------------------------------------------------------- 
    struct Dr : public std::unary_function<Edge_handle,Face_handle> { 
	Face_handle operator()(const Edge_handle& e) const { return e->dr(); }
    };
    // ------------------------------------------------------------------------- 
    struct Ul : public std::unary_function<Edge_handle,Face_handle> { 
	Face_handle operator()(const Edge_handle& e) const { return e->ul(); }
    };
    // ------------------------------------------------------------------------- 
    struct Ur : public std::unary_function<Edge_handle,Face_handle> {
	Face_handle operator()(const Edge_handle& e) const { return e->ur(); }
    };
    // -------------------------------------------------------------------------
    struct Is_xx_right 
	: public std::unary_function<Vertex_handle,bool> {
	bool operator()(const Vertex_handle& v) const { return v->is_xx_right(); }
    };
    // ------------------------------------------------------------------------- 
    struct Is_right_xx
	: public std::unary_function<Vertex_handle,bool> {
	bool operator()(const Vertex_handle& v) const { return v->is_right_xx(); }
    };
    // ------------------------------------------------------------------------- 
    struct Is_xx_left 
	: public std::unary_function<Vertex_handle,bool> {
	bool operator()(const Vertex_handle& v) const { return v->is_xx_left(); }
    };
    // ------------------------------------------------------------------------- 
    struct Is_left_xx
	: public std::unary_function<Vertex_handle,bool> {
	bool operator()(const Vertex_handle& v) const { return v->is_left_xx(); }
    };
    // ------------------------------------------------------------------------- 
};

// ----------------------------------------------------------------------------- 

template < class _Vc >
struct Visibility_complex_right_ccw_traits;

template < class _Vc >
struct Visibility_complex_left_ccw_traits;

template < class _Vc >
struct Visibility_complex_right_cw_traits;

template < class _Vc >
struct Visibility_complex_left_cw_traits;

// ----------------------------------------------------------------------------- 

template < class _Vc > 
struct Visibility_complex_ccw_traits
{
    typedef typename _Vc::Gt             Gt;
    typedef typename _Vc::Disk_handle Disk_handle;
    typedef typename _Vc::Vertex         Vertex;
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    typedef typename _Vc::Edge           Edge;
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Face_handle    Face_handle;
    typedef Visibility_complex_left_ccw_traits<_Vc>   Left_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Right_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Sup                 Sup;
    typedef typename Base::Set_sup             Set_sup;
    typedef typename Base::Inf                 Inf;
    typedef typename Base::Set_inf             Set_inf;
    typedef typename Base::Ccw_edge            Ccw_edge;
    typedef typename Base::Set_ccw_edge        Set_ccw_edge;
    typedef typename Base::Cw_edge             Cw_edge;
    typedef typename Base::Set_cw_edge         Set_cw_edge;
    typedef typename Base::Cusp_edge           Cusp_edge;
    typedef typename Base::Set_cusp_edge       Set_cusp_edge;
    typedef typename Base::Cusp_face           Cusp_face;
    // ------------------------------------------------------------------------- 
    struct Splice {
	void operator()(Edge_handle e , Vertex_handle v) const {
	    if (e == 0 || v == 0 || (e->sup() != 0 && *v == *e->sup()) || 
				    (e->inf() != 0 && *v == *e->inf())) return;
	    // splicing geometric Arc_2
	    Edge_handle tmp = new Edge(e->sign(),e->object());
	    e->split(*tmp,(e->object() == v->target_object()) ? v->target() : 
								v->source());
	    Edge_handle ep = v->ccw_edge(e->object());
	    if (ep == 0) ep = tmp; 
	    else delete tmp;
	    // Updating Edge <---> Vertex pointers 
	    ep->set_inf(v); v->set_ccw_edge(ep);
	    if (e->sup() != 0) {
		ep->set_sup(e->sup()); e->sup()->set_cw_edge(ep);
	    }
	    e->set_sup(v);   v->set_cw_edge(e);
	}
    };
    // ------------------------------------------------------------------------- 
    class Merge {
    public:
	Merge(_Vc* vc) : _vc(vc) { }
	void operator()(Edge_handle e, Edge_handle f) const {
	    // Special cases when the arc cannot be merged.
	    if (e == 0 || f == 0) return;
	    if (e->sup() == 0 || e->sup() != f->inf()) return;
	    if (e->sup()->is_constraint() || _vc->is_on_convex_hull(e->sup())) 
		return;
	    e->join(*f);
	    if (f->sup() != 0) {
		e->set_sup(f->sup());
		f->sup()->set_cw_edge(e);
	    }
	    f->inf()->set_cw_edge (0,e->object());
	    f->inf()->set_ccw_edge(0,e->object());
	    delete f;
	}
    private:
	_Vc* _vc;
    };
    // -------------------------------------------------------------------------
    /*
    struct Vertex_creator {
	Vertex_handle operator()(const Edge_handle& b , 
				 const Edge_handle& t) const {
	    typename _Vc::Type_util type;
	    return new Vertex(type(b->sign(),t->sign()),
			      b->object(),t->object());
	}
    };
    */
    // ------------------------------------------------------------------------- 
    struct Less_vertex_handle {
	bool operator()(const Vertex_handle& v1,const Vertex_handle& v2) const {
	    if (v1 == v2) return false;
	    return Less_bitangent<Gt>()(*v1,*v2);
	}
    };
    // ------------------------------------------------------------------------- 
    struct Set_adjacent_faces {
	void operator()(Edge_handle e,
			const Face_handle& f0,
			const Face_handle& f1,
			const Face_handle& f2) const { 
	    e->set_adjacent_faces(f0,f1,f2); 
	}
    };
    // ------------------------------------------------------------------------- 
    struct Set_adjacent_faces_one_to_one {
	void operator()(Edge_handle e,
			Face_handle f0, Face_handle f1, Face_handle f2) const { 
	    if (f0 != 0) f0->set_top_edge(e);
	    if (f2 != 0) f2->set_bottom_edge(e);
	    if (f1 != 0 && e->sign()) f1->set_bottom_edge(e);
	    else if (f1 != 0) f1->set_top_edge(e);
	    e->set_adjacent_faces(f0,f1,f2); 
	}
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc > 
struct Visibility_complex_cw_traits
{
    typedef typename _Vc::Gt             Gt;
    typedef typename _Vc::Disk_handle Disk_handle;
    typedef typename _Vc::Vertex         Vertex;
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    typedef typename _Vc::Edge           Edge;
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Face_handle    Face_handle;
    typedef Visibility_complex_left_cw_traits<_Vc>   Left_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>  Right_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Inf                 Sup;
    typedef typename Base::Set_inf             Set_sup;
    typedef typename Base::Sup                 Inf;
    typedef typename Base::Set_sup             Set_inf;
    typedef typename Base::Cw_edge             Ccw_edge;
    typedef typename Base::Set_cw_edge         Set_ccw_edge;
    typedef typename Base::Ccw_edge            Cw_edge;
    typedef typename Base::Set_ccw_edge        Set_cw_edge;
    // -------------------------------------------------------------------------
    struct Cusp_edge {
	Edge_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->pi()->cusp_edge(p); }
    };
    // -------------------------------------------------------------------------
    struct Set_cusp_edge {
	void operator()(Vertex_handle v, Edge_handle e ,
			const Disk_handle& p) const
	{ v->pi()->set_cusp_edge(e,p); }
    };
    // -------------------------------------------------------------------------
    struct Cusp_face {
	Face_handle operator()(const Vertex_handle& v, 
			       const Disk_handle& p) const
	{ return v->pi()->cusp_face(p); }
    };
    // -------------------------------------------------------------------------
    struct Splice {
	void operator()(Edge_handle e , Vertex_handle v) const {
	    if (e == 0 || v == 0 || (e->sup() != 0 && *v == *e->sup()) || 
				    (e->inf() != 0 && *v == *e->inf())) return;
	    // splicing geometric Arc_2
	    Edge_handle tmp = new Edge(e->sign(),e->object());
	    e->split_cw(*tmp,(e->object() == v->target_object()) ? v->target() : 
								   v->source());
	    Edge_handle ep = v->cw_edge(e->object());
	    if (ep == 0) ep = tmp; 
	    else delete tmp;
	    // Updating Edge <---> Vertex pointers 
	    ep->set_sup(v); v->set_cw_edge(ep);
	    if (e->inf() != 0) {
		ep->set_inf(e->inf()); e->inf()->set_ccw_edge(ep);
	    }
	    e->set_inf(v);   v->set_ccw_edge(e);
	}
    };
    // ------------------------------------------------------------------------- 
    class Merge {
    public:
	Merge(_Vc* vc) : _vc(vc) { }
	void operator()(Edge_handle e , Edge_handle f) const {
	    // Special cases when the arc cannot be merged.
	    if (e == 0 || f == 0) return;
	    if (e->inf() == 0 || e->inf() != f->sup()) return;
	    if (e->inf()->is_constraint() || _vc->is_on_convex_hull(e->inf())) 
		return;
	    f->join(*e); *e = *f;
	    if (f->inf() != 0) {
		e->set_inf(f->inf());
		f->inf()->set_ccw_edge(e);
	    }
	    f->sup()->set_ccw_edge (0,e->object());
	    f->sup()->set_cw_edge(0,e->object());
	    delete f;
	}
    private:
	_Vc* _vc;
    };
    // -------------------------------------------------------------------------
    /*
    struct Vertex_creator {
	Vertex_handle operator()(const Edge_handle& b , 
				 const Edge_handle& t) const {
	    typename _Vc::Type_util type;
	    return new Vertex(type(t->sign(),b->sign()),
			      t->object(),b->object());
	}
    };
    */
    // ------------------------------------------------------------------------- 
    struct Less_vertex_handle {
	bool operator()(const Vertex_handle& v1,const Vertex_handle& v2) const{
	    if (v1 == v2) return false;
	    return Greater_bitangent<Gt>()(*v1,*v2);
	}
    };
    // ------------------------------------------------------------------------- 
    struct Set_adjacent_faces {
	void operator()(Edge_handle e,
			const Face_handle& f0,
			const Face_handle& f1,
			const Face_handle& f2) const { 
	    if (e->sign()) e->set_adjacent_faces(f0,f2,f1); 
	    else e->set_adjacent_faces(f1,f0,f2);
	}
    };
    // ------------------------------------------------------------------------- 
    struct Set_adjacent_faces_one_to_one {
	void operator()(Edge_handle e,
			Face_handle f0,Face_handle f1,Face_handle f2) const { 
	    if (e->sign()) {
		if (f0 != 0) f0->set_top_edge(e);
		if (f1 != 0) f1->set_bottom_edge(e);
		if (f2 != 0) f2->set_bottom_edge(e);
		e->set_adjacent_faces(f0,f2,f1); 
	    }
	    else {
		if (f0 != 0) f0->set_top_edge(e);
		if (f1 != 0) f1->set_top_edge(e);
		if (f2 != 0) f2->set_bottom_edge(e);
		e->set_adjacent_faces(f1,f0,f2); 
	    }
	}
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc >
struct Visibility_complex_source_target_traits
{
    // -------------------------------------------------------------------------
    typedef typename _Vc::Vertex_handle Vertex_handle;
    typedef typename _Vc::Vertex        Vertex;
    typedef typename _Vc::Edge_handle   Edge_handle;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Source_object Source_object;
    typedef typename Base::Target_object Target_object;
    // -------------------------------------------------------------------------
    struct Vertex_creator {
	Vertex_handle operator()(const Edge_handle& b , 
				 const Edge_handle& t) const {
	    typename _Vc::Type_util type;
	    return new Vertex(type(b->sign(),t->sign()),
			      b->object(),t->object());
	}
    };
    // -------------------------------------------------------------------------
};

template < class _Vc >
struct Visibility_complex_target_source_traits
{
    // -------------------------------------------------------------------------
    typedef typename _Vc::Vertex_handle Vertex_handle;
    typedef typename _Vc::Vertex        Vertex;
    typedef typename _Vc::Edge_handle   Edge_handle;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Target_object Source_object;
    typedef typename Base::Source_object Target_object;
    // -------------------------------------------------------------------------
    struct Vertex_creator {
	Vertex_handle operator()(const Edge_handle& b , 
				 const Edge_handle& t) const {
	    typename _Vc::Type_util type;
	    return new Vertex(type(t->sign(),b->sign()),
			      t->object(),b->object());
	}
    };
    // -------------------------------------------------------------------------
};

template < class _Vc >
struct Visibility_complex_right_traits
{
    // -------------------------------------------------------------------------
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Top_edge    Top_edge;
    typedef typename Base::Bottom_edge Bottom_edge;
    // -------------------------------------------------------------------------
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Face_handle    Face_handle;
    // -------------------------------------------------------------------------
    struct Back_view
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& f) const 
	    { return f->back_view(); }
    };
    // -------------------------------------------------------------------------
    struct Front_view
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& f) const 
	    { return f->front_view(); }
    };
    // -------------------------------------------------------------------------
    struct Sign
	: public std::unary_function<Edge_handle,bool> {
	bool operator()(const Edge_handle& e) const { return e->sign(); }
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc >
struct Visibility_complex_left_traits
{
    // -------------------------------------------------------------------------
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Bottom_edge  Top_edge;
    typedef typename Base::Top_edge     Bottom_edge;
    // -------------------------------------------------------------------------
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Face_handle    Face_handle;
    // -------------------------------------------------------------------------
    struct Back_view
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& f) const 
	    { return f->front_view(); }
    };
    // -------------------------------------------------------------------------
    struct Front_view
	: public std::unary_function<Face_handle,Edge_handle> {
	Edge_handle operator()(const Face_handle& f) const 
	    { return f->back_view(); }
    };
    // ------------------------------------------------------------------------- 
    struct Sign
	: public std::unary_function<Edge_handle,bool> {
	bool operator()(const Edge_handle& e) const { return !e->sign(); }
    };
    // ------------------------------------------------------------------------- 
};

template < class C_ , class S_ , class R_ >
struct Visibility_complex_traits : public C_ , public S_ , public R_
{
    typedef typename C_::Disk_handle Disk_handle;
    typedef typename C_::Vertex_handle  Vertex_handle;
    typedef typename C_::Vertex         Vertex;
    typedef typename C_::Edge           Edge;
    typedef typename C_::Edge_handle    Edge_handle;
    typedef typename C_::Face_handle    Face_handle;
    typedef typename C_::Ccw_edge       Ccw_edge;
    typedef typename C_::Cw_edge        Cw_edge;
    typedef typename C_::Cusp_edge      Cusp_edge;
    typedef typename C_::Set_cusp_edge  Set_cusp_edge;
    typedef typename C_::Cusp_face      Cusp_face;
    typedef typename C_::Sup            Sup;
    typedef typename C_::Inf            Inf;
    typedef typename C_::Less_vertex_handle Less_vertex_handle;
    //typedef typename C_::Vertex_creator Vertex_creator;
    typedef typename S_::Target_object  Target_object;
    typedef typename S_::Source_object  Source_object;
    typedef typename S_::Vertex_creator Vertex_creator;
    typedef typename R_::Sign           Sign;
    // -------------------------------------------------------------------------
    struct Ccw_target_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Ccw_edge()(v,Target_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Ccw_source_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Ccw_edge()(v,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Cw_target_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Cw_edge()(v,Target_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Cw_source_edge
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Cw_edge()(v,Source_object()(v)); }
    };
    // ------------------------------------------------------------------------- 
    struct Set_target_cusp_edge {
	void operator()(Vertex_handle v,Edge_handle e) const
	{ Set_cusp_edge()(v,e,Target_object()(v)); }
    };
    // ------------------------------------------------------------------------- 
    struct Set_source_cusp_edge {
	void operator()(Vertex_handle v,Edge_handle e) const
	{ Set_cusp_edge()(v,e,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Target_cusp_edge 
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Cusp_edge()(v,Target_object()(v)); }
    };
    // ------------------------------------------------------------------------- 
    struct Source_cusp_edge 
	: public std::unary_function<Vertex_handle,Edge_handle> {
	Edge_handle operator()(const Vertex_handle& v) const
	{ return Cusp_edge()(v,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Target_cusp_face 
	: public std::unary_function<Vertex_handle,Face_handle> {
	Face_handle operator()(const Vertex_handle& v) const
	{ return Cusp_face()(v,Target_object()(v)); }
    };
    // ------------------------------------------------------------------------- 
    struct Source_cusp_face 
	: public std::unary_function<Vertex_handle,Face_handle> {
	Face_handle operator()(const Vertex_handle& v) const
	{ return Cusp_face()(v,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Cc
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v,
				 const Disk_handle& p) const { 
	    typename C_::Ccw_edge ccw_edge; typename C_::Sup sup;
	    return (ccw_edge(v,p) == 0) ?  0 : sup(ccw_edge(v,p)); 
	}
    };
    // -------------------------------------------------------------------------
    struct Cw
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v,
				 const Disk_handle& p) const { 
	    typename C_::Cw_edge cw_edge; typename C_::Inf inf;
	    return (cw_edge(v,p) == 0) ? 0 : inf(cw_edge(v,p)); 
	}
    };
    // -------------------------------------------------------------------------
    struct CcR
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v) const
	{ return Cc()(v,Target_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct CcL
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v) const
	{ return Cc()(v,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct CwR
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v) const
	{ return Cw()(v,Source_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct CwL
	: public std::unary_function<Vertex_handle,Vertex_handle> {
	Vertex_handle operator()(const Vertex_handle& v) const
	{ return Cw()(v,Target_object()(v)); }
    };
    // -------------------------------------------------------------------------
    struct Walk_vertex {
	Vertex_handle operator()(const Edge_handle& b , 
				 const Edge_handle& t) const {
	    Sup sup; Inf inf; Vertex_creator bit; Sign sign; 
	    Vertex_handle tmp = bit(b,t); Less_vertex_handle chi2;
	    if ((sup(b) != 0 && (!sign(b) || sup(b)->is_constraint()) && 
				 chi2(sup(b),tmp)) ||
		(sup(t) != 0 && (sign(t)  || sup(t)->is_constraint()) &&
				 chi2(sup(t),tmp))) { delete tmp; return 0; }
	    if (sup(b) != 0 && *tmp == *sup(b)) { delete tmp; return sup(b); }
	    if (inf(b) != 0 && *tmp == *inf(b)) { delete tmp; return inf(b); }
	    if (sup(t) != 0 && *tmp == *sup(t)) { delete tmp; return sup(t); }
	    if (inf(t) != 0 && *tmp == *inf(t)) { delete tmp; return inf(t); }
	    return tmp;
	}
    };
    // -------------------------------------------------------------------------
};

template < class _Vc >
struct Visibility_complex_right_ccw_traits
: public Visibility_complex_traits<Visibility_complex_ccw_traits<_Vc>,
				   Visibility_complex_source_target_traits<_Vc>,
				   Visibility_complex_right_traits<_Vc> >
{
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_left_ccw_traits<_Vc>   Left_ccw_traits;
    typedef Visibility_complex_left_ccw_traits<_Vc>   Left_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Right_ccw_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Right_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Left_cw_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Right_cw_traits;
    typedef Right_cw_traits                           Cw_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Dl Dl;
    typedef typename Base::Dr Dr;
    typedef typename Base::Ul Ul;
    typedef typename Base::Ur Ur;
    typedef typename Base::Is_xx_left Is_xx_left;
    typedef typename Base::Is_left_xx Is_left_xx;
    // ------------------------------------------------------------------------- 
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    struct Chi3 {
	bool operator()(const Vertex_handle& v, const Vertex_handle& w) const {
	    typename _Vc::Gt::Orientation_infinite chi3;
	    return (chi3(*v,*w) == LEFTTURN);
	}
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc >
struct Visibility_complex_left_ccw_traits
: public Visibility_complex_traits<Visibility_complex_ccw_traits<_Vc> ,
				   Visibility_complex_target_source_traits<_Vc>,
				   Visibility_complex_left_traits<_Vc> >
{
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_left_ccw_traits<_Vc>   Right_ccw_traits;
    typedef Visibility_complex_left_ccw_traits<_Vc>   Right_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Left_ccw_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Left_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Right_cw_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Left_cw_traits;
    typedef Left_cw_traits                            Cw_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Ur Dl;
    typedef typename Base::Ul Dr;
    typedef typename Base::Dr Ul;
    typedef typename Base::Dl Ur;
    typedef typename Base::Is_right_xx Is_xx_left;
    typedef typename Base::Is_xx_right Is_left_xx;
    // ------------------------------------------------------------------------- 
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    struct Chi3 {
	bool operator()(const Vertex_handle& v, const Vertex_handle& w) const {
	    typename _Vc::Gt::Orientation_infinite chi3;
	    return (chi3(*v,*w->pi()) == RIGHTTURN);
	}
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc >
struct Visibility_complex_right_cw_traits
: public Visibility_complex_traits<Visibility_complex_cw_traits<_Vc> ,
				   Visibility_complex_target_source_traits<_Vc>,
				   Visibility_complex_right_traits<_Vc> >
{
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_left_ccw_traits<_Vc>   Left_cw_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Right_cw_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Left_ccw_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Left_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Right_ccw_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Right_traits;
    typedef Right_cw_traits                           Cw_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Dr Dl;
    typedef typename Base::Dl Dr;
    typedef typename Base::Ur Ul;
    typedef typename Base::Ul Ur;
    typedef typename Base::Is_left_xx Is_xx_left;
    typedef typename Base::Is_xx_left Is_left_xx;
    // ------------------------------------------------------------------------- 
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    struct Chi3 {
	bool operator()(const Vertex_handle& v, const Vertex_handle& w) const {
	    typename _Vc::Gt::Orientation_infinite chi3;
	    return (chi3(*v,*w->pi()) == LEFTTURN);
	}
    };
    // ------------------------------------------------------------------------- 
};

template < class _Vc >
struct Visibility_complex_left_cw_traits
: public Visibility_complex_traits<Visibility_complex_cw_traits<_Vc> ,
				   Visibility_complex_source_target_traits<_Vc>,
				   Visibility_complex_left_traits<_Vc> >
{
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_left_ccw_traits<_Vc>   Right_cw_traits;
    typedef Visibility_complex_right_ccw_traits<_Vc>  Left_cw_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Right_ccw_traits;
    typedef Visibility_complex_left_cw_traits<_Vc>    Right_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Left_ccw_traits;
    typedef Visibility_complex_right_cw_traits<_Vc>   Left_traits;
    typedef Left_cw_traits                            Cw_traits;
    // ------------------------------------------------------------------------- 
    typedef Visibility_complex_duality_traits_base<_Vc>   Base;
    typedef typename Base::Ul Dl;
    typedef typename Base::Ur Dr;
    typedef typename Base::Dl Ul;
    typedef typename Base::Dr Ur;
    typedef typename Base::Is_xx_right Is_xx_left;
    typedef typename Base::Is_right_xx Is_left_xx;
    // ------------------------------------------------------------------------- 
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    struct Chi3 {
	bool operator()(const Vertex_handle& v, const Vertex_handle& w) const {
	    typename _Vc::Gt::Orientation_infinite chi3;
	    return (chi3(*v,*w) == RIGHTTURN);
	}
    };
    // ------------------------------------------------------------------------- 
};

// ----------------------------------------------------------------------------- 

// Functions to convert traits
template <class _Vc>
Visibility_complex_right_cw_traits<_Vc> 
change_ccw_cw(const Visibility_complex_right_ccw_traits<_Vc>&)
{ return Visibility_complex_right_cw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_right_ccw_traits<_Vc> 
change_ccw_cw(Visibility_complex_right_cw_traits<_Vc>) 
{ return Visibility_complex_right_ccw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_left_cw_traits<_Vc> 
change_ccw_cw(Visibility_complex_left_ccw_traits<_Vc>)
{ return Visibility_complex_left_cw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_left_ccw_traits<_Vc> 
change_ccw_cw(Visibility_complex_left_cw_traits<_Vc>)
{ return Visibility_complex_left_ccw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_right_cw_traits<_Vc> 
change_right_left(Visibility_complex_left_cw_traits<_Vc>)  
{ return Visibility_complex_right_cw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_right_ccw_traits<_Vc> 
change_right_left(Visibility_complex_left_ccw_traits<_Vc>)
{ return Visibility_complex_right_ccw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_left_cw_traits<_Vc> 
change_right_left(Visibility_complex_right_cw_traits<_Vc>)
{ return Visibility_complex_left_cw_traits<_Vc>(); }

template <class _Vc>
Visibility_complex_left_ccw_traits<_Vc> 
change_right_left(Visibility_complex_right_ccw_traits<_Vc>)
{ return Visibility_complex_left_ccw_traits<_Vc>(); }


// ----------------------------------------------------------------------------- 

CGAL_END_NAMESPACE

#endif //  VISIBILITY_COMPLEX_CCW_CW_TRAITS_H
