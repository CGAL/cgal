#ifndef VISIBILITY_COMPLEX_SWEEP_ITERATOR_H
#define VISIBILITY_COMPLEX_SWEEP_ITERATOR_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <list>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr >
class Visibility_complex_sweep_iterator_base
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_sweep_iterator_base<_Vc,_Tp,_Ref,_Ptr> Self;
    // -------------------------------------------------------------------------
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t            difference_type;
    typedef _Tp   value_type;
    typedef _Ptr  pointer;
    typedef _Ref  reference;
    // -------------------------------------------------------------------------
    Visibility_complex_sweep_iterator_base() : a(0) , min(0) { }
    Visibility_complex_sweep_iterator_base(_Vc* ant , pointer m) 
	: a(ant) , min(m) { }
    Visibility_complex_sweep_iterator_base(const Self& vi) 
	: a(vi.a), min(vi.min) { }
    Visibility_complex_sweep_iterator_base& operator=(const Self& vi)
    {
	if (this != &vi) {
	    a   = vi.a;
	    min = vi.min;
	}
	return *this;
    }
    // -------------------------------------------------------------------------
    bool operator==(const Self& vi) const {return (a == vi.a && min == vi.min);}
    bool operator!=(const Self& vi) const {return !(*this == vi);}
    // -------------------------------------------------------------------------
    reference operator*()  const { return *min; }
    pointer   operator->() const { return (min == 0) ? 0 : &(*min); }
    // -------------------------------------------------------------------------
protected:
    _Vc* a;
    pointer min;
};

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr >
class Visibility_complex_sweep_iterator
    : public Visibility_complex_sweep_iterator_base<_Vc,_Tp,_Ref,_Ptr>
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_sweep_iterator_base<_Vc,_Tp,_Ref,_Ptr> Base;
    typedef Visibility_complex_sweep_iterator<_Vc,_Tp,_Ref,_Ptr> Self;

  typedef typename Base::pointer pointer;
    // -------------------------------------------------------------------------
    Visibility_complex_sweep_iterator() : Base() { }
    Visibility_complex_sweep_iterator(_Vc* ant , pointer m) : Base(ant,m) { }
    Visibility_complex_sweep_iterator(_Vc* ant) : Base(ant,0) 
    { 
	min = non_swept_minimal();
    }
    Visibility_complex_sweep_iterator(const Self& vi) : Base(vi) { }
    // -------------------------------------------------------------------------
    Self operator++() { 
	if (min == 0) return *this;
	a->sweep(min);
	min = non_swept_minimal();
	return *this; 
    }
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    typename Base::pointer non_swept_minimal() {
	typename _Vc::Minimals_iterator m = a->minimals_begin();
	while (m != a->minimals_end() && a->is_swept(&(*m))) 
	    a->erase_minimal(m++);
	return (m != a->minimals_end())? &(*m) : 0;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc , class _Tp , class _Ref , class _Ptr , class _Is_upward >
class Visibility_complex_linear_sweep_iterator
    : public Visibility_complex_sweep_iterator_base<_Vc,_Tp,_Ref,_Ptr>
{
public:
    // -------------------------------------------------------------------------
    typedef Visibility_complex_sweep_iterator_base<_Vc,_Tp,_Ref,_Ptr> Base;
    typedef Visibility_complex_linear_sweep_iterator<_Vc,_Tp,_Ref,_Ptr,
						     _Is_upward>       Self;
    typedef typename Base::pointer pointer;
    // -------------------------------------------------------------------------
    Visibility_complex_linear_sweep_iterator() : Base() { }
    Visibility_complex_linear_sweep_iterator(_Vc* ant , pointer m) 
	: Base(ant,m) { }
    Visibility_complex_linear_sweep_iterator(_Vc* ant) : Base(ant,0) 
    { 
	min = non_swept_minimal();
    }
    Visibility_complex_linear_sweep_iterator(const Self& vi) : Base(vi) { }
    // -------------------------------------------------------------------------
    Self operator++() { 
	if (min == 0) return *this;
	a->sweep(min);
	min = non_swept_minimal();
	return *this; 
    }
    // -------------------------------------------------------------------------
private:
    // -------------------------------------------------------------------------
    pointer non_swept_minimal() {
	_Is_upward is_upward;
	typename _Vc::Minimals_iterator m = a->minimals_begin();
	while (m != a->minimals_end() && !is_upward(*m))
	    a->erase_minimal(m++);
	return (m != a->minimals_end())? &(*m) : 0;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_SWEEP_ITERATOR_H
