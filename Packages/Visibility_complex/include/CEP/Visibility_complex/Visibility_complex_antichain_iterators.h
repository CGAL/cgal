#ifndef VISIBILITY_COMPLEX_ANTICHAIN_ITERATORS_H
#define VISIBILITY_COMPLEX_ANTICHAIN_ITERATORS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class _Ant >
class Vc_antichain_iterator_base
{
public:
    // -------------------------------------------------------------------------
    typedef typename _Ant::Edge_const_iterator  Edge_const_iterator;
    typedef typename _Ant::Face_handle          Face_handle;
    typedef Vc_antichain_iterator_base<_Ant>    Self;
    // -------------------------------------------------------------------------
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t            difference_type;
    // -------------------------------------------------------------------------
    Vc_antichain_iterator_base() : _f(0), _a(0) , _e(0) { }
    Vc_antichain_iterator_base(const _Ant* a , Edge_const_iterator e) 
	: _a(a) , _e(e)
    {
	if (_e != _a->edges_end()) {
	    _f = _e->dl();
	    if (_f == 0 || _f == _a->infinite_face()) increment();
	}
	else _f = _a->infinite_face();
    }
    Vc_antichain_iterator_base(const Self& fi)
	: _f(fi._f), _a(fi._a) , _e(fi._e) { }
    Vc_antichain_iterator_base& operator=(const Self& fi)
    {
	_f       = fi._f;
	_a       = fi._a;
	_e       = fi._e;
	return *this;
    }
    // -------------------------------------------------------------------------
    bool operator==(const Self& fi) const {return (_a == fi._a && _e == fi._e);}
    bool operator!=(const Self& fi) const {return !(*this == fi);}
    // -------------------------------------------------------------------------
    void increment()
    {
	if (_e == _a->edges_end()) return;
	Face_handle old = _f;
	_f = _e->dr();
	if (_f == 0 || _f == _a->infinite_face() || _f == old) 
	{ ++_e; _f = _e->dl(); }
	if (_f == 0 || _f == _a->infinite_face() || _f == old) increment();
    }
    Self operator++() { increment(); return *this; }
    // -------------------------------------------------------------------------
protected:
    Face_handle   _f;
    const _Ant*   _a;
    Edge_const_iterator _e;
};

// -----------------------------------------------------------------------------

template < class _Ant , class _Tp , class _Ref , class _Ptr >
class Vc_antichain_face_iterator 
    : public Vc_antichain_iterator_base<_Ant>
{
public:
    // -------------------------------------------------------------------------
    typedef Vc_antichain_iterator_base<_Ant>                Base;
    typedef Vc_antichain_face_iterator<_Ant,_Tp,_Ref,_Ptr>  Self;    
    typedef typename Base::Edge_const_iterator  Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Vc_antichain_face_iterator() : Base() { }
    Vc_antichain_face_iterator(const _Ant* a , Edge_const_iterator e) 
	: Base(a,e) { }
    // -------------------------------------------------------------------------
    reference operator*()  const { return *_f; }
    pointer   operator->() const { return _f;  }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Ant , class _Tp , class _Ref , class _Ptr , class _Sup >
class Vc_antichain_vertex_iterator 
    : public Vc_antichain_iterator_base<_Ant>
{
public:
    // -------------------------------------------------------------------------
    typedef Vc_antichain_iterator_base<_Ant>  Base;
    typedef typename Base::Edge_const_iterator               Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Vc_antichain_vertex_iterator() : Base() { }
    Vc_antichain_vertex_iterator(const _Ant* a , Edge_const_iterator e)
	: Base(a,e) { }
    // -------------------------------------------------------------------------
    _Ref operator*() { 
	_Sup sup;
#ifdef ITERATOR_DEBUG
	cout << "_f = " << long(_f) << " sup(_f) = " << long(sup(_f)) << endl;
#endif
	while (sup(_f) == 0 && _e != _a->edges_end()) increment();
	return *sup(_f);
    }
    _Ptr operator->() { 
	_Sup sup;
#ifdef ITERATOR_DEBUG
	cout << "_f = " << long(_f) << " sup(_f) = " << long(sup(_f)) << endl;
#endif
	while (sup(_f) == 0 && _e != _a->edges_end()) increment();
	return sup(_f);
    }
    // -------------------------------------------------------------------------
    Vc_antichain_vertex_iterator operator++() { 
	//_Sup sup;
	//while (sup(_f) == 0 && _e != _a->edges_end()) increment(); 
	//if (_e != _a->edges_end()) _f = 0;
	increment();
	return *this; 
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
