// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_SWEEP_ITERATOR_H
#define CGAL_VISIBILITY_COMPLEX_2_SWEEP_ITERATOR_H

#include <CGAL/basic.h>

#include <list>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template < class Vc_ , class Tp_ , class Ref_ , class Ptr_ >
class Sweep_iterator_base
{
public:
    // -------------------------------------------------------------------------
    typedef Sweep_iterator_base<Vc_,Tp_,Ref_,Ptr_> Self;
    // -------------------------------------------------------------------------
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t            difference_type;
    typedef Tp_   value_type;
    typedef Ptr_  pointer;
    typedef Ref_  reference;
    // -------------------------------------------------------------------------
    Sweep_iterator_base() : a(0) , min(0) { }
    Sweep_iterator_base(Vc_* ant , pointer m) 
	: a(ant) , min(m) { }
    Sweep_iterator_base(const Self& vi) 
	: a(vi.a), min(vi.min) { }
    Sweep_iterator_base& operator=(const Self& vi)
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
    Vc_* a;
    pointer min;
};

// -----------------------------------------------------------------------------

template < class Vc_ , class Tp_ , class Ref_ , class Ptr_ >
class Sweep_iterator
    : public Sweep_iterator_base<Vc_,Tp_,Ref_,Ptr_>
{
public:
    // -------------------------------------------------------------------------
    typedef Sweep_iterator_base<Vc_,Tp_,Ref_,Ptr_> Base;
    typedef Sweep_iterator<Vc_,Tp_,Ref_,Ptr_> Self;

  typedef typename Base::pointer pointer;

  using Base::min;
  using Base::a;
    // -------------------------------------------------------------------------
    Sweep_iterator() : Base() { }
    Sweep_iterator(Vc_* ant , pointer m) : Base(ant,m) { }
    Sweep_iterator(Vc_* ant) : Base(ant,0) 
    { 
	min = non_swept_minimal();
    }
    Sweep_iterator(const Self& vi) : Base(vi) { }
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
	typename Vc_::Minimals_iterator m = a->minimals_begin();
	while (m != a->minimals_end() && m->final_antichain()// a->is_swept(&(*m))
               ) 
	    a->erase_minimal(m++);
	return (m != a->minimals_end())? &(*m) : 0;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class Vc_ , class Tp_ , class Ref_ , class Ptr_ , class Is_upward_ >
class Linear_sweep_iterator
    : public Sweep_iterator_base<Vc_,Tp_,Ref_,Ptr_>
{
public:
    // -------------------------------------------------------------------------
    typedef Sweep_iterator_base<Vc_,Tp_,Ref_,Ptr_> Base;
    typedef Linear_sweep_iterator<Vc_,Tp_,Ref_,Ptr_,
						     Is_upward_>       Self;
    typedef typename Base::pointer pointer;
    // -------------------------------------------------------------------------
  using Base::min;
  using Base::a;

    Linear_sweep_iterator() : Base() { }
    Linear_sweep_iterator(Vc_* ant , pointer m) 
	: Base(ant,m) { }
    Linear_sweep_iterator(Vc_* ant) : Base(ant,0) 
    { 
	min = non_swept_minimal();
    }
    Linear_sweep_iterator(const Self& vi) : Base(vi) { }
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
	Is_upward_ is_upward;
	typename Vc_::Minimals_iterator m = a->minimals_begin();
	while (m != a->minimals_end() && !is_upward(*m))
	    a->erase_minimal(m++);
	return (m != a->minimals_end())? &(*m) : 0;
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_SWEEP_ITERATOR_H
