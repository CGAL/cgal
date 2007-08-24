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

#ifndef CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_ITERATORS_H
#define CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_ITERATORS_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------

template < class Ant_ >
class Antichain_iterator_base
{
public:
    // -------------------------------------------------------------------------
    typedef typename Ant_::Edge_const_iterator  Edge_const_iterator;
    typedef typename Ant_::Face_const_handle          Face_const_handle;
    typedef Antichain_iterator_base<Ant_>    Self;
    // -------------------------------------------------------------------------
    typedef std::forward_iterator_tag iterator_category;
    typedef std::ptrdiff_t            difference_type;
    // -------------------------------------------------------------------------
    Antichain_iterator_base() : f_(0), a_(0) , e_(0) { }
    Antichain_iterator_base(const Ant_* a , Edge_const_iterator e) 
	: a_(a) , e_(e)
    {
	if (e_ != a_->edges_end()) {
	    f_ = e_->dl();
	    if (f_ == 0 || f_ == a_->infinite_face()) increment();
	}
	else f_ = a_->infinite_face();
    }
    Antichain_iterator_base(const Self& fi)
	: f_(fi.f_), a_(fi.a_) , e_(fi.e_) { }
    Antichain_iterator_base& operator=(const Self& fi)
    {
	f_       = fi.f_;
	a_       = fi.a_;
	e_       = fi.e_;
	return *this;
    }
    // -------------------------------------------------------------------------
    bool operator==(const Self& fi) const {return (a_ == fi.a_ && e_ == fi.e_);}
    bool operator!=(const Self& fi) const {return !(*this == fi);}
    // -------------------------------------------------------------------------
    void increment()
    {
	if (e_ == a_->edges_end()) return;
	Face_const_handle old = f_;
	f_ = e_->dr();
	if (f_ == 0 || f_ == a_->infinite_face() || f_ == old) 
	{ ++e_; f_ = e_->dl(); }
	if (f_ == 0 || f_ == a_->infinite_face() || f_ == old) increment();
    }
    Self operator++() { increment(); return *this; }
    // -------------------------------------------------------------------------
protected:
    Face_const_handle   f_;
    const Ant_*   a_;
    Edge_const_iterator e_;
};

// -----------------------------------------------------------------------------

template < class Ant_ , class Tp_ , class Ref_ , class Ptr_ >
class Antichain_face_iterator 
    : public Antichain_iterator_base<Ant_>
{
public:
    // -------------------------------------------------------------------------
    typedef Antichain_iterator_base<Ant_>                Base;
    typedef Antichain_face_iterator<Ant_,Tp_,Ref_,Ptr_>  Self;    
    typedef typename Base::Edge_const_iterator  Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef Tp_   value_type;
    typedef Ptr_  pointer;
    typedef Ref_  reference;
    // -------------------------------------------------------------------------
    Antichain_face_iterator() : Base() { }
    Antichain_face_iterator(const Ant_* a , Edge_const_iterator e) 
	: Base(a,e) { }
    // -------------------------------------------------------------------------
    reference operator*()  const { return const_cast<reference>(*this->f_); }
    pointer   operator->() const { return const_cast<pointer>(this->f_);  }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class Ant_ , class Tp_ , class Ref_ , class Ptr_ , class Sup_ >
class Antichain_vertex_iterator 
    : public Antichain_iterator_base<Ant_>
{
public:
    // -------------------------------------------------------------------------
    typedef Antichain_iterator_base<Ant_>  Base;
    typedef typename Base::Edge_const_iterator               Edge_const_iterator;
    // -------------------------------------------------------------------------
    typedef Tp_   value_type;
    typedef Ptr_  pointer;
    typedef Ref_  reference;

    using Base::increment;
    // -------------------------------------------------------------------------
    Antichain_vertex_iterator() : Base() { }
    Antichain_vertex_iterator(const Ant_* a , Edge_const_iterator e)
	: Base(a,e) { }
    // -------------------------------------------------------------------------
    Ref_ operator*() { 
	Sup_ sup;
	while (sup(this->f_) == 0 && this->e_ != this->a_->edges_end())
          increment();
	return const_cast<reference>(*sup(this->f_));
    }
    Ptr_ operator->() { 
	Sup_ sup;
	while (sup(this->f_) == 0 && this->e_ != this->a_->edges_end())
          increment();
	return const_cast<pointer>(sup(this->f_));
    }
    // -------------------------------------------------------------------------
    Antichain_vertex_iterator operator++() { 
	//Sup_ sup;
	//while (sup(f_) == 0 && e_ != a_->edges_end()) increment(); 
	//if (e_ != a_->edges_end()) f_ = 0;
	increment();
	return *this; 
    }
};

// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_ITERATORS_H
