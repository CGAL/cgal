// Copyright (c) 2006-2011  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_MODIFIABLE_PRIORITY_QUEUE_H
#define CGAL_MODIFIABLE_PRIORITY_QUEUE_H

#include <climits> // Needed by the following Boost header for CHAR_BIT.
#include <boost/optional.hpp>

#ifdef CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
#include <boost/pending/relaxed_heap.hpp>
#else
#include <CGAL/STL_Extension/internal/boost/mutable_queue.hpp>
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
  #else
  typedef  internal::boost_::mutable_queue<IndexedType,std::vector<IndexedType>,Compare,ID> Heap;
  #endif //CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;

public:

  Modifiable_priority_queue( size_type largest_ID, Compare const& c = Compare(), ID const& id = ID() ) : mHeap(largest_ID,c,id) {}

  void push ( value_type const& v ) { mHeap.push(v) ; }

  void update ( value_type const& v ) { mHeap.update(v); }

  void erase ( value_type const& v  ) { mHeap.remove(v); }

  value_type top() const { return mHeap.top(); }

  void pop() { mHeap.pop(); }

  bool empty() const { return mHeap.empty() ; }

  bool contains ( value_type const& v ) { return mHeap.contains(v) ; }

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

  value_type top_and_pop()
  {
    CGAL_precondition(!empty());
    value_type v = top();
    pop();
    return v;
  }

private:

  Heap mHeap ;

} ;

} //namespace CGAL

#endif

