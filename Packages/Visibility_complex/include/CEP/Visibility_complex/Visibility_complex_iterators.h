#ifndef VISIBILITY_COMPLEX_ITERATORS_H
#define VISIBILITY_COMPLEX_ITERATORS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <list>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template<class _Vc>
class Visibility_complex_iterator_base
{
public:
    // -------------------------------------------------------------------------
    typedef typename _Vc::Vertex		  Vertex;
    typedef typename _Vc::Vertex_handle		  Vertex_handle;
    typedef typename std::list<Vertex_handle>::iterator    iterator;
    typedef typename std::list<Vertex_handle>::const_iterator const_iterator;
    typedef Visibility_complex_iterator_base<_Vc> Self;
    // -------------------------------------------------------------------------
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t            difference_type;
    // -------------------------------------------------------------------------
    Visibility_complex_iterator_base() : _v(0) {}
    Visibility_complex_iterator_base(iterator v) : _v(v) {}
    Visibility_complex_iterator_base(const Self& vi) : _v(vi._v) {}
    Visibility_complex_iterator_base& operator=(const Self& vi)
    {
	_v       = vi._v;
	return *this;
    }
    // -------------------------------------------------------------------------
    bool operator==(const Self& vi) const { return (_v == vi._v);  }
    bool operator!=(const Self& vi) const { return !(*this == vi); }
    // -------------------------------------------------------------------------
    Vertex&       operator*()  const { return *(*_v);    }
    Vertex_handle operator->() const { return &(*(*_v)); }
    // -------------------------------------------------------------------------
    void increment() { ++_v; }
    // -------------------------------------------------------------------------
protected:
    iterator _v;
};

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr >
class Visibility_complex_vertex_iterator
    : public Visibility_complex_iterator_base<_Vc> 
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_iterator_base<_Vc>        Base;
    typedef typename Base::iterator                               iterator;
    typedef Visibility_complex_vertex_iterator<_Vc,_Tp,_Ref,_Ptr> Self;
    // -------------------------------------------------------------------------
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Visibility_complex_vertex_iterator() : Base()  {}
    Visibility_complex_vertex_iterator(iterator v) : Base(v) {}
    Visibility_complex_vertex_iterator(Base v) : Base(v) {}
    // -------------------------------------------------------------------------
    Self operator++() { increment(); return *this; }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr >
class Visibility_complex_edge_iterator
    : public Visibility_complex_iterator_base<_Vc>
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_iterator_base<_Vc>        Base;
    typedef typename Base::iterator                               iterator;
    typedef Visibility_complex_edge_iterator<_Vc,_Tp,_Ref,_Ptr> Self;
    // -------------------------------------------------------------------------
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Visibility_complex_edge_iterator() : Base(),visited(false)  {}
    Visibility_complex_edge_iterator(iterator v) : Base(v) {}
    Visibility_complex_edge_iterator(Base v) : Base(v) ,visited(false) {}
    Visibility_complex_edge_iterator& operator=(const Self& vi)
    {
	Base::operator=(vi);
	visited  = vi.visited;
	return *this;
    }
    // -------------------------------------------------------------------------
    bool operator==(const Self& vi) const 
	{ return (Base::operator==(vi) && visited == vi.visited);  }
    // -------------------------------------------------------------------------
    Self operator++() { 
	if (visited) { increment(); visited = false; }
	else visited = true;
	return *this;
    }
    // -------------------------------------------------------------------------
    reference operator*()  const
	{ return (visited) ? *(Base::operator->()->ccw_target_edge()) : 
			     *(Base::operator->()->ccw_source_edge());}
    pointer   operator->() const
	{ return (visited) ? Base::operator->()->ccw_target_edge()    : 
			     Base::operator->()->ccw_source_edge(); }
    // -------------------------------------------------------------------------
private:
    bool visited;
};

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr >
struct Visibility_complex_face_iterator
    : public Visibility_complex_iterator_base<_Vc>
{
    // -------------------------------------------------------------------------
    typedef Visibility_complex_iterator_base<_Vc>        Base;
    typedef typename Base::iterator                               iterator;
    typedef Visibility_complex_face_iterator<_Vc,_Tp,_Ref,_Ptr> Self;
    // -------------------------------------------------------------------------
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Visibility_complex_face_iterator()       : Base()  {} 
    Visibility_complex_face_iterator(iterator v) : Base(v) {}
    Visibility_complex_face_iterator(Base v) : Base(v) {} 
    // -------------------------------------------------------------------------
    Self operator++() { increment(); return *this; }
    // -------------------------------------------------------------------------
    reference operator*()  const { return *(Base::operator->()->inf()); }
    pointer   operator->() const { return Base::operator->()->inf();    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
