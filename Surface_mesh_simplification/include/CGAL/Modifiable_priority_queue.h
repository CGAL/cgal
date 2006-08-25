// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_MODIFIABLE_PRIORITY_QUEUE_H
#define CGAL_MODIFIABLE_PRIORITY_QUEUE_H

#include <boost/relaxed_heap.hpp>


CGAL_BEGIN_NAMESPACE


template <class Type_ 
         ,class Compare_ = CGAL::Compare<Type_>
         ,class Allocator_ = CGAL_ALLOCATOR(int)
         >
class Modifiable_priority_queue
{
  
  template<class T, class Compare>
  struct Less_from_compare
  {
    Less_from_compare( Compare const& c ) : mComp(c) {}
    
    bool operator() ( T const& a, T const& b ) const
    {
      return mComp(a,b) == SMALLER ;
    }
    
    Compare mComp ;
  } ;
  
  template<class T>
  struct Id_map : public boost::put_get_helper<std::size_t, Id_map<T> >
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t                      value_type;
    typedef std::size_t                      reference;
    typedef T                                key_type;
  
    Id_map() {}
  
    reference operator[](key_type const& e) const 
    {
      return e->ID ;
    }
  };
  
public:

  
  typedef boost::relaxed_heap<Type_,Less_from_compare<Type_,Compare_>,Id_map<Type_> > Heap ;
  
  typedef Modifiable_priority_queue Self;
  
  typedef Compare_ Compare;
  
  typedef typename Heap::value_type      value_type;
  typedef typename Heap::size_type       size_type;
  
  struct handle 
  {
    handle() {}
    handle ( value_type const& v_ ) : v(v_) {} 
    
    friend bool operator == ( handle const& x, handle const& y ) { return x.v == y.v ; }
    friend bool operator != ( handle const& x, handle const& y ) { return !(x==y); }
    value_type v ; 
  } ;
  
public:

  Modifiable_priority_queue( size_type n, Compare const& c ) : mHeap(n, Less_from_compare<Type_,Compare_>(c) ) {}
  
  handle push ( value_type const& v ) { mHeap.push(v) ; return handle(v) ; }
  
  handle update ( handle h ) { mHeap.update(h.v); return h ; }
  
  void erase ( handle h ) { mHeap.remove(h.v); }

  value_type top() const { return mHeap.top() ; }
  
  void pop() { mHeap.pop(); }
  
  bool empty() const { return mHeap.empty() ; }
  
  size_type size() const { return mHeap.size() ; }

  //void set_compare ( Compare const& c ) { mHeap.key_comp() = c ; }
  
  boost::optional<value_type> extract_top()
  {
    boost::optional<value_type> r ;
    if ( !empty() )
    {
      value_type v = top();
      pop();
      r = boost::optional<value_type>(v) ;
    }  
    return r ;
  }
  
private:

  Heap mHeap ;  
    
} ;

CGAL_END_NAMESPACE

#endif
 
