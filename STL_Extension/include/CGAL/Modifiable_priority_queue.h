// Copyright (c) 2006-2011  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_MODIFIABLE_PRIORITY_QUEUE_H
#define CGAL_MODIFIABLE_PRIORITY_QUEUE_H

#include <climits> // Neeeded by the following Boost header for CHAR_BIT.
#include <boost/optional.hpp>
#ifdef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
#include <boost/pending/relaxed_heap.hpp>
#else
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/graph/properties.hpp>
#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP

namespace CGAL {

template <class IndexedType_ 
         ,class Compare_ = std::less<IndexedType_>
         ,class ID_      = boost::identity_property_map
         >
class Modifiable_priority_queue
{
public:

  typedef Modifiable_priority_queue Self;
  
  typedef IndexedType_ IndexedType ;
  typedef Compare_     Compare;
  typedef ID_          ID ;
  
#ifdef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  typedef boost::relaxed_heap<IndexedType,Compare,ID> Heap;
  typedef bool handle ;

  Modifiable_priority_queue( size_type largest_ID, Compare const& c, ID const& id ) : mHeap(largest_ID,c,id) {}
  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;
  handle push ( value_type const& v ) { mHeap.push(v) ; return handle(true) ; }
  handle erase ( value_type const& v, handle  ) { mHeap.remove(v); return null_handle() ; }
  handle erase ( value_type const& v  ) { mHeap.remove(v); return null_handle() ; }
  bool contains ( value_type const& v ) { return mHeap.contains(v) ; }
  handle update ( value_type const& v, handle h ) { mHeap.update(v); return h ; }
  static handle null_handle() { return handle(false); }
#else
  struct Reverse_compare{
    const Compare c;
    Reverse_compare(){}
    Reverse_compare(Compare const& c):c(c){}
    template<typename T>
     bool operator() (T const& a, T const& b) const
     {
       return !c(a,b);
     }
  };
  typedef boost::heap::fibonacci_heap<IndexedType,boost::heap::compare<Reverse_compare> > Heap;
  typedef typename Heap::handle_type handle;
  //the fibonacci_heap uses the inverse of the compare used in the relaxed heap.
  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;
public:
  Modifiable_priority_queue( size_type, Compare const& c, ID const& ):mHeap(Reverse_compare(c)) {}
  handle push ( value_type const& v ) { return mHeap.push(v) ;}
  handle erase ( value_type const& v, handle h) { mHeap.erase(h); return null_handle() ; }
  handle update ( value_type const& v, handle h ) { mHeap.update(h); return h ; }
  static handle null_handle() { return handle(); }
#endif //CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP


  value_type top() const { return mHeap.top() ; }
  
  void pop() { mHeap.pop(); }
  
  bool empty() const { return mHeap.empty() ; }

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

} //namespace CGAL

#endif
 
