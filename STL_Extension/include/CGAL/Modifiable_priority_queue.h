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

#include <CGAL/STL_Extension/internal/boost/relaxed_heap.hpp>
#include <CGAL/STL_Extension/internal/boost/mutable_queue.hpp>

#include <boost/heap/fibonacci_heap.hpp>

#include <type_traits>

namespace CGAL {

enum Heap_type { CGAL_BOOST_FIBONACCI_HEAP, CGAL_BOOST_PENDING_MUTABLE_QUEUE, CGAL_BOOST_PENDING_RELAXED_HEAP };

template <class IndexedType_
         ,class Compare_ = std::less<IndexedType_>
         ,class ID_      = boost::identity_property_map
         , Heap_type heap_type = CGAL_BOOST_PENDING_MUTABLE_QUEUE
         >
class Modifiable_priority_queue
{
public:

  typedef Modifiable_priority_queue Self;

  typedef IndexedType_ IndexedType ;
  typedef Compare_     Compare;
  typedef ID_          ID ;

  typedef std::conditional_t<heap_type == CGAL_BOOST_PENDING_MUTABLE_QUEUE,
                             internal::boost_::mutable_queue<IndexedType,std::vector<IndexedType>,Compare,ID>,
                             internal::boost_::relaxed_heap<IndexedType,Compare,ID> > Heap;

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

template <class IndexedType_
         ,class Compare_ 
         ,class ID_>
class Modifiable_priority_queue<IndexedType_, Compare_, ID_, CGAL_BOOST_FIBONACCI_HEAP>
{
public:

  typedef Modifiable_priority_queue Self;

  typedef IndexedType_ IndexedType ;
  typedef Compare_     Compare;
  typedef ID_          ID ;

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

  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;

public:

  Modifiable_priority_queue( size_type largest_ID, Compare const& c = Compare(), ID const& id = ID() )
    : mHeap(Reverse_compare(c))
    , mID(id)
    , mHandles(largest_ID)
  {}

  void push ( value_type const& v ) { mHandles[get(mID, v)]=mHeap.push(v) ; }

  void update ( value_type const& v ) { mHeap.update(mHandles[get(mID, v)]); }

  void erase ( value_type const& v  ) { 
    auto vid = get(mID, v);
    mHeap.erase(mHandles[vid]);
    mHandles[vid]=typename Heap::handle_type();
  }

  value_type top() const { return mHeap.top(); }

  void pop() { mHeap.pop(); }

  bool empty() const { return mHeap.empty() ; }

  bool contains ( value_type const& v ) { return mHandles[get(mID, v)] != typename Heap::handle_type(); }

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
  ID mID;
  std::vector<typename Heap::handle_type> mHandles;
} ;

} //namespace CGAL

#endif

