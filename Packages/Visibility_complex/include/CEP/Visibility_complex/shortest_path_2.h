#ifndef SHORTEST_PATH_2_H
#define SHORTEST_PATH_2_H

#ifndef VISIBILITY_COMPLEX_2_H
#include <CEP/Visibility_complex/Visibility_complex_2.h>
#endif

#include <functional>
#include <set>
#include <math.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class _Vc >
class Sh_atom 
{
public:
    typedef typename _Vc::Vertex_handle  Vertex_handle;
    typedef typename _Vc::Edge_handle    Edge_handle;
    typedef typename _Vc::Gt::Exact_NT   Exact_NT;
    typedef typename _Vc::Disk_handle Disk_handle;
    typedef Sh_atom<_Vc>                 Atom;
private:
    Vertex_handle _v;
    Edge_handle   _e;
public:
    Sh_atom() : _v(0) , _e(0) {}
    Sh_atom(Vertex_handle v) : _v(v) , _e(0) {}
    Sh_atom(Edge_handle e) : _v(0) , _e(e) {}
    bool operator==(const Atom& a) const
	{ return (_e == a._e && _v == a._v); }
    bool operator!=(const Atom& a) const
	{ return !(*this == a); }

    Edge_handle   edge()        const { return _e; }
    Vertex_handle vertex()      const { return _v; }

    long index()                const 
	{ return (_e == 0) ? long(_v)          : long(_e);       }
    Exact_NT distance()         const 
	{ return (_e == 0) ? _v->distance()    : _e->distance(); }
    Exact_NT weight()           const 
	{ return (_e == 0) ? _v->weight()      : _e->weight();   }
    Atom   prev()               const 
	{ return (_e == 0) ? _v->prev()        : _e->prev() ;    }
    Vertex_handle next_vertex(Disk_handle s) const 
	{ return (_e == 0) ? _v->next_vertex(s) : _e->next_vertex(s); }
    Edge_handle next_edge(Disk_handle s)     const 
	{ return (_e == 0) ? _v->next_edge(s)   : _e->next_edge(s); }
};

// -----------------------------------------------------------------------------

template < class _Vc >
class Sh_edge
    : public Visibility_complex_edge_base<_Vc>
{
public:
    typedef typename _Vc::Gt                           Gt;
    typedef typename _Vc::Disk_handle               Disk_handle;
    typedef typename Gt::Exact_NT                      Exact_NT;
    typedef typename Gt::Arc_2                  Arc;
    typedef typename Gt::Point_2                       Point_2;
    typedef Visibility_complex_edge_base<_Vc>          Base;
    typedef typename Base::Vertex_handle                        Vertex_handle;
    typedef typename Base::Edge_handle                          Edge_handle;
    //typedef typename Arc::Vertex_iterator       Arc_iterator;
    //typedef typename Arc::Vertex_const_iterator Arc_const_iterator;
    typedef Sh_atom<_Vc>                               Atom;

private:
    Exact_NT      _distance;
    Atom          _prev;

public:
    // -------------------------------------------------------------------------
    Sh_edge() : Base() , _distance(Exact_NT(-1)) { }
    Sh_edge(bool s,Disk_handle p) : Base(s,p), _distance(Exact_NT(-1)) { }
    Sh_edge(Vertex_handle v0 , Vertex_handle v1 , Disk_handle p)
	: Base(v0,v1,p) , _distance(Exact_NT(-1)) { }
    // -------------------------------------------------------------------------
    Exact_NT weight() const { return Gt().length(*this,*inf(),*sup()); }
    // -------------------------------------------------------------------------
    Exact_NT distance() const { return _distance; }
    void   set_distance(Exact_NT d) { _distance = d; }
    // -------------------------------------------------------------------------
    Atom  prev() const { return _prev; }
    void  set_prev(Atom v) { _prev = v; }
    // -------------------------------------------------------------------------
    Vertex_handle next_vertex(Disk_handle s) const
    {
	if (object() != s) {
	    if ((sign()  && sup()->target_object() == object()) ||
	        (!sign() && inf()->target_object() == object())) return 0;
	}
	return (sign()) ? sup() : inf();
    }
    Edge_handle   next_edge(Disk_handle s)   const
    {
	if (object() != s) {
	    if (( sign() && sup()->is_constraint()) || 
		(!sign() && inf()->is_constraint())) return 0;
	}
	if (sign()) return (sup()->target_object() == object()) ?
			sup()->ccw_target_edge()  : sup()->ccw_source_edge();
	return (inf()->target_object() == object()) ?
	    inf()->cw_target_edge() : inf()->cw_source_edge();
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc >
class Sh_vertex
    : public Visibility_complex_vertex_base<_Vc>
{
public:
    typedef typename _Vc::Gt                       Gt;
    typedef typename _Vc::Gt::Exact_NT             Exact_NT;
    typedef typename _Vc::Vertex_handle            Vertex_handle;
    typedef typename _Vc::Edge_handle              Edge_handle;
    typedef typename _Vc::Bitangent_2              Bitangent_2;
    typedef Visibility_complex_vertex_base<_Vc>    Base;
    typedef typename Base::Disk_handle                   Disk_handle;
    typedef typename Base::Type                             Type;
    typedef Sh_atom<_Vc>                           Atom;

private:
    Exact_NT _distance;
    Atom     _prev;

public:
    // -------------------------------------------------------------------------
    Sh_vertex() : Base() , _distance(Exact_NT(-1)) { }
    Sh_vertex(Type t , Disk_handle start , Disk_handle finish)
	: Base (t,start,finish) , _distance(Exact_NT(-1)) { }
    Sh_vertex(Edge_handle start , Edge_handle finish)
	: Base(start,finish)    , _distance(Exact_NT(-1)) { }
    Sh_vertex(const Bitangent_2& b) : Base(b) , _distance(Exact_NT(-1)) { }
    // -------------------------------------------------------------------------
    Exact_NT weight()   const { return Gt().length(*this); }
    // -------------------------------------------------------------------------
    Exact_NT distance() const { return _distance; }
    void   set_distance(Exact_NT d) { _distance = d; }
    // -------------------------------------------------------------------------
    Atom prev() const { return _prev; }
    void set_prev(Atom e) { _prev = e; }
    // -------------------------------------------------------------------------
  Vertex_handle next_vertex(Disk_handle /*s*/) const { return 0; }
    Edge_handle   next_edge(Disk_handle s)   const {
	if (is_constraint() && s != source_object() && s != target_object()) {
	    if (is_left_right() && prev().edge() ==  cw_source_edge()) return 0;
	    if (is_right_left() && prev().edge() == ccw_source_edge()) return 0;
	}
	return (is_xx_left()) ? ccw_target_edge() : cw_target_edge(); 
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

class Sh_items : public Visibility_complex_items {
public:
    template <class _Vc >
    struct Edge_wrapper {
	typedef Sh_edge<_Vc>   Edge;
    };
    template <class _Vc>
    struct Vertex_wrapper {
	typedef Sh_vertex<_Vc> Vertex;
    };
};

// -----------------------------------------------------------------------------

template <class At>
struct Less_atom {
    bool operator() (const At& a, const At& b) const {
	typedef typename At::Exact_NT Exact_NT;
	if (a.distance() == b.distance()) return (a.index() < b.index());
	if (a.distance() == Exact_NT(-1)) return false;
	if (b.distance() == Exact_NT(-1)) return true;
	return (a.distance() < b.distance());
    }
};

// -----------------------------------------------------------------------------
// Computes all the shortest paths from Edge s

template < class _Edge_handle , class _Gtr , class _It>
void
__shortest_path_2(_Edge_handle s , Visibility_complex_2<_Gtr,_It>* /*V*/) // warning V never used
{
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Exact_NT                     Exact_NT;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_antichain<_Gtr,_It>      Antichain;
    typedef typename Antichain::Vertex                           Vertex;
    typedef Visibility_complex_2<_Gtr,_It>              Visibility_complex;
    typedef typename Visibility_complex::Vertex_iterator         Vertex_iterator;
    typedef Sh_atom<Antichain>                          Atom;
    // -------------------------------------------------------------------------
    Atom start(s);
    // The priority queue, we push start ---------------------------------------
    typedef std::set<Atom, Less_atom<Atom> > Queue;
    Queue X;    
    start.edge()->set_distance(0);
    X.insert(start);
    // Dijkstra algorithm ------------------------------------------------------
    while (!X.empty()) 
    {
	Atom a = *X.begin(); 
	X.erase(X.begin());

	Exact_NT da = a.distance(); 
	typename Antichain::Vertex_handle v = a.next_vertex(s->object());
	typename Antichain::Edge_handle   e = a.next_edge(s->object());

	if (v != 0) {
	    Exact_NT d = da + v->weight();
	    if (v->distance() == -1 || d < v->distance()) {
		if (v->distance() != Exact_NT(-1)) {
		    typename Queue::iterator vit = X.find(Atom(v));
		    if (vit != X.end()) X.erase(vit);
		}
		v->set_distance(d); v->set_prev(a);
		X.insert(Atom(v));
	    }
	}
	if (e != 0) {
	    Exact_NT d = da + e->weight();
	    if (e->distance() == -1 || d < e->distance()) {
		if (e->distance() != Exact_NT(-1)) {
		    typename Queue::iterator eit = X.find(Atom(e));
		    if (eit != X.end()) X.erase(eit);
		}
		e->set_distance(d); e->set_prev(a);
		X.insert(Atom(e));
	    }
	} 
    }
}

// -----------------------------------------------------------------------------
// Assumes __shortest_path_2(s) has been called.

template < class _Edge_handle , class OutputIterator , class _Gtr , class _It>
OutputIterator
__recover_path(_Edge_handle s , _Edge_handle t , 
	       OutputIterator result , Visibility_complex_2<_Gtr,_It>* /*V*/) // warning: V never used
{
    typedef Visibility_complex_antichain<_Gtr,_It>  Antichain;
    typedef typename Antichain::Vertex                       Vertex;
    typedef Sh_atom<Antichain>                      Atom;

    Atom start(s); Atom finish(t);
    std::list<Vertex> path;
    while ((finish.vertex() != 0 || finish.edge() != 0) && finish != start) {
	if (finish.vertex() != 0) path.push_back(*finish.vertex());
	finish = finish.prev();
    }
    if (finish == start) copy(path.begin(),path.end(),result);
    return result;
}

// -----------------------------------------------------------------------------

template < class _Edge_handle , class OutputIterator , class _Gtr , class _It>
OutputIterator
__shortest_path_2(_Edge_handle s , _Edge_handle t , 
		  OutputIterator result , Visibility_complex_2<_Gtr,_It>* V)
{
    __shortest_path_2(s,V);
    __recover_path(s,t,result,V);
    return result;
}

// -----------------------------------------------------------------------------
// Compute the shortest path from the Edge s to all the edges on object t

template < class OutputIterator, class _Gtr , class _It>
typename _Gtr::Exact_NT
__recover_path(typename Visibility_complex_2<_Gtr,_It>::Edge_handle s,
	       typename Visibility_complex_2<_Gtr,_It>::Disk_handle t,
	       OutputIterator result, Visibility_complex_2<_Gtr,_It>* V)
{
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Exact_NT                     Exact_NT;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_antichain<_Gtr,_It>      Antichain;
    typedef typename Antichain::Vertex                           Vertex;
    typedef Visibility_complex_2<_Gtr,_It>              Visibility_complex;
    typedef typename Antichain::Edge_iterator                    Edge_iterator;
    typedef typename Antichain::Edge_handle                      Edge_handle;
    // Stop if two disks intersect ---------------------------------------------
    if (!V->is_valid()) return Exact_NT(-1);
    // Find edges two edges with opposite sign on t ----------------------------
    typename Antichain::Edge_handle ept(0),emt(0);
    Edge_iterator e = V->antichain()->edges_begin();
    while (ept == 0 || emt == 0) {
	CGAL_precondition(e != V->antichain()->edges_end());
	if (e->object() == t &&  e->sign() && ept == 0) ept = &(*e);
	if (e->object() == t && !e->sign() && emt == 0) emt = &(*e);
	++e;
    }
    // closest edge to s on t+  ------------------------------------------------
    Edge_handle min = ept;
    Edge_handle f   = ept; 
    do {
	if (f->distance() < min->distance()) min = f;
	f = f->sup()->ccw_edge(f->object());
    } while (f != ept);
    // closest edge to s on t- -------------------------------------------------
    f = emt; 
    do {
	if (f->distance() < min->distance()) min = f;
	f = f->sup()->ccw_edge(f->object());
    } while (f != emt);
    // Shortest path from s to object t ----------------------------------------
    __recover_path(s,min,result,V);
    return min->distance();
}

//------------------------------------------------------------------------------

template < class OutputIterator, class _Gtr , class _It>
OutputIterator
shortest_path_2(typename Visibility_complex_2<_Gtr,_It>::Disk_handle s,
		typename Visibility_complex_2<_Gtr,_It>::Disk_handle t,
		OutputIterator result, Visibility_complex_2<_Gtr,_It>* V)
{
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Exact_NT                     Exact_NT;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_antichain<_Gtr,_It>      Antichain;
    typedef typename Antichain::Vertex                           Vertex;
    typedef Visibility_complex_2<_Gtr,_It>              Visibility_complex;
    typedef typename Antichain::Edge_iterator                    Edge_iterator;
    typedef typename Antichain::Edge_handle                      Edge_handle;
    // Stop is V if two disks intersect ----------------------------------------
    if (!V->is_valid()) return result;
    // Find two edge on s with opposite sign -----------------------------------
    typename Antichain::Edge_handle eps(0),ems(0);
    Edge_iterator e = V->antichain()->edges_begin();
    while (eps == 0 || ems == 0) {
	CGAL_precondition(e != V->antichain()->edges_end());
	if (e->object() == s &&  e->sign() && eps == 0) eps = &(*e);
	if (e->object() == s && !e->sign() && ems == 0) ems = &(*e);
	++e;
    }
    // Shortest path from eps to t ---------------------------------------------
    __shortest_path_2(eps,V);
    std::list<Vertex> pathp;
    Exact_NT minp_ft = __recover_path(eps,t,back_inserter(pathp),V);
    // Shortest path from ems to t ---------------------------------------------
    __shortest_path_2(ems,V);
    std::list<Vertex> pathm;
    Exact_NT minm_ft = __recover_path(ems,t,back_inserter(pathm),V);
    // Compare the two paths and return the smallest ---------------------------
    if (minm_ft < minp_ft) copy(pathm.begin(),pathm.end(),result);
    else                   copy(pathp.begin(),pathp.end(),result);
    return result;
}

// -----------------------------------------------------------------------------

template < class InputIterator , class OutputIterator, class _Gtr , class _It>
OutputIterator
all_shortest_path_2(InputIterator first, InputIterator last,
		    typename Visibility_complex_2<_Gtr,_It>::Disk_handle s,
		    OutputIterator result, Visibility_complex_2<_Gtr,_It>* V)
{
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Exact_NT                 Exact_NT;
    typedef Visibility_complex_antichain<_Gtr,_It>  Antichain;
    typedef typename Antichain::Vertex                       Vertex;
    typedef typename Antichain::Disk_handle               Disk_handle;
    // Stop if two disks intersect ---------------------------------------------
    if (!V->is_valid()) return result;
    // Find edges on s ---------------------------------------------------------
    typename Antichain::Edge_handle ep(0),em(0);
    typename Antichain::Edge_iterator ei = V->antichain()->edges_begin();
    while (ep == 0 || em == 0) {
	CGAL_precondition(ei != V->antichain()->edges_end());
	if (ei->object() == s &&  ei->sign() && ep == 0) ep = &(*ei);
	if (ei->object() == s && !ei->sign() && em == 0) em = &(*ei);
	++ei;
    }
    CGAL_precondition(ep != 0 && em != 0);
    // Finding all the points of the scene -------------------------------------
    std::vector<Disk_handle> points;
    for (InputIterator d = first; d != last ; ++d)
	if (&(*d) != s && _Gtr().is_vertex(*d)) points.push_back(&(*d));
    // Shortest path from ep to all points -------------------------------------
    __shortest_path_2(ep,V);
    std::vector<std::list<Vertex> > pathp(points.size());
    std::vector<Exact_NT>           pft;
    int i = 0;
    for (typename std::vector<Disk_handle>::iterator eh  = points.begin(); 
					       eh != points.end() ; ++eh , ++i) 
	pft.push_back(__recover_path(ep,*eh,back_inserter(pathp[i]),V));
    // Shortest path from em to all points -------------------------------------
    __shortest_path_2(em,V);
    std::vector<std::list<Vertex> > pathm(points.size());
    std::vector<Exact_NT>           mft;
    i = 0;
    for (typename std::vector<Disk_handle>::iterator eh  = points.begin(); 
					       eh != points.end() ; ++eh , ++i) 
	mft.push_back(__recover_path(em,*eh,back_inserter(pathm[i]),V));
    // Returning shortest path to all points -----------------------------------
    for (unsigned int i = 0 ; i < points.size() ; i++) {
	if (mft[i] < pft[i]) copy(pathm[i].begin(),pathm[i].end(),result);
	else                 copy(pathp[i].begin(),pathp[i].end(),result);
    }
    return result;
}

// -----------------------------------------------------------------------------

template < class InputIterator , class ConstraintIt, class OutputIterator , 
	   class _Gtr , class _It>
OutputIterator
shortest_path_2(InputIterator first, InputIterator last ,
		ConstraintIt  firstc,ConstraintIt  lastc,
		const typename _Gtr::Point_2& p, 
		const typename _Gtr::Point_2& q,
		OutputIterator result, _Gtr , _It )
{
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Exact_NT                     Exact_NT;
    typedef typename _Gtr::Point_2                      Point_2;
    typedef typename _Gtr::Disk                         Disk;
    typedef Disk*                                       Disk_handle;
    // -------------------------------------------------------------------------
    typedef Visibility_complex_antichain<_Gtr,_It>      Antichain;
    typedef typename Antichain::Vertex                           Vertex;
    typedef Visibility_complex_2<_Gtr,_It>              Visibility_complex;
    // -------------------------------------------------------------------------
    
    // Add the points p and q to the scene -------------------------------------
    std::list<Disk> objects;
    std::map<long,Disk_handle> Disk_map;
    for (InputIterator it = first; it != last; it++) {
	objects.push_back(*it);
	Disk_map[long(&(*it))] = &objects.back();
    }
    std::list<Vertex> L;
    for (ConstraintIt c = firstc; c != lastc; c++)
	L.push_back(Vertex(c->type(),Disk_map[long(c->source_object())],
				     Disk_map[long(c->target_object())]));
    objects.push_back(_Gtr().make_convex_from_point(p)); 
    Disk_handle s = &objects.back();
    objects.push_back(_Gtr().make_convex_from_point(q)); 
    Disk_handle t = &objects.back();
    // -------------------------------------------------------------------------
    // Compute visibility complex.
    Visibility_complex* V = new Visibility_complex(objects.begin(),
						   objects.end(),
						   L.begin(),L.end());
    // -------------------------------------------------------------------------
    shortest_path_2(s,t,result,V);
    return result;
}

// -----------------------------------------------------------------------------

template < class InputIterator , class OutputIterator , 
	   class _Gtr , class _It>
OutputIterator
shortest_path_2(InputIterator first, InputIterator last ,
		const typename _Gtr::Point_2& p, 
		const typename _Gtr::Point_2& q,
		OutputIterator result, _Gtr , _It )
{
    typedef Visibility_complex_2<_Gtr,_It> Vc;
    std::list<typename Vc::Vertex> V;
    shortest_path_2(first,last,V.begin(),V.end(), p,q,result,
		    _Gtr(),_It());  // does not get resolved
    return result;
}

// -----------------------------------------------------------------------------

template < class InputIterator , class ConstraintIt, 
	   class OutputIterator, class _Gtr >
OutputIterator
shortest_path_2(InputIterator first, InputIterator last ,
		ConstraintIt  firstc,ConstraintIt  lastc,
		const typename _Gtr::Point_2& p, 
		const typename _Gtr::Point_2& q,
		OutputIterator result, _Gtr )
{
    shortest_path_2(first,last,firstc,lastc, p,q,result,
		    _Gtr(),Sh_items());
    return result;
}

// -----------------------------------------------------------------------------

template < class InputIterator , class OutputIterator , class _Gtr >
OutputIterator
shortest_path_2(InputIterator first, InputIterator last ,
		const typename _Gtr::Point_2& p, 
		const typename _Gtr::Point_2& q,
		OutputIterator result, _Gtr )
{
    shortest_path_2(first,last,p,q,result, _Gtr(),Sh_items());
    return result;
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
