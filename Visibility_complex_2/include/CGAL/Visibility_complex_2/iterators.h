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

#ifndef CGAL_VISIBILITY_COMPLEX_2_ITERATORS_H
#define CGAL_VISIBILITY_COMPLEX_2_ITERATORS_H

#include <CGAL/basic.h>

#include <vector>
#include <iterator>

CGAL_BEGIN_NAMESPACE

namespace Visibility_complex_2_details {

struct Which_edge_wrapper {
  enum Which {CCW_SOURCE,CCW_TARGET,SOURCE_CUSP,TARGET_CUSP};
};

struct Which_face_wrapper {
  enum Which {UNSET,INF,SOURCE_CUSP,TARGET_CUSP};
};

template<class VertexBasicIterator> struct Bogus_iterator {
  struct Base {
    VertexBasicIterator v_;
    bool opposite;
  };
  struct Vertex {
    typedef typename Bogus_iterator::Base Base;
  };
  struct Edge :Base,Which_edge_wrapper {
    typedef typename Bogus_iterator::Base Base;
    Which which;
  };
  struct Face :Base, Which_face_wrapper {
    typedef typename Bogus_iterator::Base Base;
    Which which;
  };
};



template<class Vc_,class VertexBasicIterator,class Ptr_,
         class From_=typename Bogus_iterator<VertexBasicIterator>::Base>
class Iterator_base 
{
public:
    typedef VertexBasicIterator                   iterator;
private:
    typedef Iterator_base<Vc_,VertexBasicIterator,Ptr_,From_> Self;
public:
    Iterator_base() : v_(0),opposite(true) {}
    Iterator_base(const iterator& v) : v_(v),opposite(true) {}
    Iterator_base(const Self& vi) : v_(vi.v_),opposite(true) {}
    Iterator_base(const From_& fi) : v_(fi.v_),opposite(fi.opposite) {}
    Iterator_base& operator=(const Self& vi)
    {
	v_       = vi.v_;
        opposite = vi.opposite;
	return *this;
    }

    bool operator==(const Self& vi) const { 
      return (v_ == vi.v_)&&opposite==vi.opposite;
    }
    bool operator!=(const Self& vi) const { return !(*this == vi); }

    Ptr_ operator->() const {
      if (opposite) return ((*v_)->pi()); else return (*v_); 
    }

    void increment() { 
      if (!opposite) {
        ++v_;
        opposite=true;
      } else opposite=false;
    }
    void decrement() { 
      if (opposite) {
        --v_;
        opposite=false;
      } else opposite=true;
    }
private:
    iterator v_;
    bool opposite;
    template<class a,class b,class c,class d> friend class Iterator_base;
};



template < class Base_,
           class Tp_,class Ref_,class Ptr_,
           class From_=
             typename Bogus_iterator<typename Base_::iterator>::Vertex>
class Vertex_iterator
  : private Base_,
    public std::iterator<std::bidirectional_iterator_tag,
                         Tp_,ptrdiff_t,Ptr_,Ref_>
{
    typedef Base_ Base;
    typedef typename Base::iterator                               iterator;
    typedef Vertex_iterator<Base,Tp_,Ref_,Ptr_,From_> Self;

    using Base::increment;
    using Base::decrement;

    typedef
      std::iterator<std::bidirectional_iterator_tag,Tp_,ptrdiff_t,Ptr_,Ref_>
      it;
public:
    Vertex_iterator() : Base()  {}
    Vertex_iterator(const iterator& v) : Base(v) {}
    Vertex_iterator(const Base& v) : Base(v) {}
    Vertex_iterator(const From_&v) : Base(v) {}
    Self& operator++() { increment(); return *this; }
    Self& operator--() { decrement(); return *this; }
    Self operator++(int) {
      Self s(*this);
      ++*this;
      return s;
    }
    Self operator--(int) {
      Self s(*this);
      --*this;
      return s;
    }
    typename it::reference operator*()  const { return *(Base::operator->()); }
    typename it::pointer   operator->() const { return Base::operator->();    }

    bool operator==(const Self& vi) const { 
      return this->Base::operator==(vi);
    }
    bool operator!=(const Self& vi) const { 
      return !(*this==vi);
    }
    template<class a,class b,class c,class d,class e>
      friend class Vertex_iterator;
};


template < class Base_,
           class Tp_,class Ref_,class Ptr_,
           class From_=
             typename Bogus_iterator<typename Base_::iterator>::Edge >
class Edge_iterator
    : private Base_,
      public std::iterator<std::bidirectional_iterator_tag,
                           Tp_,ptrdiff_t,Ptr_,Ref_>,
      private Which_edge_wrapper
{
    typedef Base_ Base;
    typedef typename Base::iterator                               iterator;
    typedef Edge_iterator<Base_,Tp_,Ref_,Ptr_,From_>  Self;

    using Base::increment;
    using Base::decrement;
    typedef
      std::iterator<std::bidirectional_iterator_tag,Tp_,ptrdiff_t,Ptr_,Ref_>
      it;

public:
    Edge_iterator() : Base(),which(CCW_SOURCE)  {}
    Edge_iterator(const iterator& v) : Base(v), which(CCW_SOURCE) {}
    Edge_iterator(const Base& v) : Base(v), which(CCW_SOURCE) {}
    Edge_iterator(const Self& v) : Base(v), which(v.which) {}
    Edge_iterator(const From_& v) : Base(v), which(v.which) {}
    Edge_iterator& operator=(const Self& vi)
    {
	Base::operator=(vi);
	which  = vi.which;
	return *this;
    }

    bool operator==(const Self& vi) const 
	{ return (Base::operator==(vi) && which == vi.which);  }
    bool operator!=(const Self& vi) const { 
      return !(*this==vi);
    }

    Self& operator++() { 
      if (which==(Base::operator->()->is_constraint()?TARGET_CUSP:CCW_TARGET)){
        increment();
        which=CCW_SOURCE;
      }	else which=static_cast<Which>(which+1);
	return *this;
    }
    Self& operator--() { 
      if (which==CCW_SOURCE){
        decrement();
        which=(Base::operator->()->is_constraint()?TARGET_CUSP:CCW_TARGET);
      }	else which=static_cast<Which>(which-1);
      return *this;
    }
    Self operator++(int) {
      Self s(*this);
      ++*this;
      return s;
    }
    Self operator--(int) {
      Self s(*this);
      --*this;
      return s;
    }
    // -------------------------------------------------------------------------
    typename it::reference operator*()  const {
      return *operator->();
    }
    typename it::pointer   operator->() const { 
      switch (which) {
      case CCW_SOURCE: return Base::operator->()->ccw_source_edge();
      case CCW_TARGET: return Base::operator->()->ccw_target_edge();
      case SOURCE_CUSP: return Base::operator->()->source_cusp_edge();
      case TARGET_CUSP: return Base::operator->()->target_cusp_edge();
      }
      std::cerr<<"which="<<which<<"\n";
      CGAL_assertion(false);
      return 0;
    }
private:
    Which which;
    template<class a,class b,class c,class d,class e>
      friend class Edge_iterator;
};


template < class Base_,
           class Tp_,class Ref_,class Ptr_,
           class From_ =
             typename Bogus_iterator<typename Base_::iterator>::Face>
class Face_iterator
    : private Base_,
      public std::iterator<std::bidirectional_iterator_tag,
                           Tp_,ptrdiff_t,Ptr_,Ref_>,
      private Which_face_wrapper
{
public:
    typedef Base_ Base;
private:
    typedef typename Base::iterator                               iterator;
    typedef Face_iterator<Base,Tp_,Ref_,Ptr_,From_>  Self;

    using Base::increment;
    using Base::decrement;
    typedef
      std::iterator<std::bidirectional_iterator_tag,Tp_,ptrdiff_t,Ptr_,Ref_>
      it;

public:
    Which which;
private:
    void set_which() const {
      (const_cast<Self*>(this))->which=
        (Base::operator->()->is_constraint()?SOURCE_CUSP:INF);
    }
public:
    Face_iterator()       : Base()  {} 
    Face_iterator(const iterator& v) : 
      Base(v),which(UNSET) {}
    Face_iterator(const Base& v) :
      Base(v),which(UNSET) {} 
    Face_iterator(const From_& v) : Base(v), which(v.which) {}
    bool operator==(const Self& fi) const {
      if (Base::operator==(fi)) {
        if (which!=UNSET) {
          if (fi.which==UNSET) fi.set_which();
        } else if (fi.which!=UNSET) set_which(); else return true;
        return which==fi.which;
      }
      return false;
    }
    bool operator!=(const Self& vi) const { 
      return !(*this==vi);
    }
    Self& operator++() {
      if (which==UNSET) set_which();
      if (which==SOURCE_CUSP)
        which=TARGET_CUSP;
      else {
        increment(); 
        which=UNSET;
      }
      return *this; }
    Self& operator--() {
      if (which==TARGET_CUSP)
        which=SOURCE_CUSP;
      else {
        decrement();
        which=(Base::operator->()->is_constraint()?TARGET_CUSP:INF);
      }
      return *this;
    }
    Self operator++(int) {
      Self s(*this);
      ++*this;
      return s;
    }
    Self operator--(int) {
      Self s(*this);
      --*this;
      return s;
    }
    typename it::reference operator*()  const {
      return *(Base::operator->()->inf());
    }
    typename it::pointer   operator->() const {
      if (which==UNSET) set_which();
      switch (which) {
      case INF: return Base::operator->()->inf();
      case SOURCE_CUSP: return Base::operator->()->source_cusp_face();
      case TARGET_CUSP: return Base::operator->()->target_cusp_face();
      default: CGAL_assertion(false);
      }
      return 0;
    }
    template<class a,class b,class c,class d,class e>
      friend class Face_iterator;
};


}
CGAL_END_NAMESPACE

#endif //CGAL_VISIBILITY_COMPLEX_2_ITERATORS_H
