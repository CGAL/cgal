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
#include <optional>

#include <CGAL/STL_Extension/internal/boost/relaxed_heap.hpp>
#include <CGAL/STL_Extension/internal/boost/mutable_queue.hpp>

#include <boost/heap/pairing_heap.hpp>

#include <algorithm>
#include <type_traits>

namespace CGAL {

enum Heap_type { CGAL_BOOST_PAIRING_HEAP, CGAL_BOOST_PENDING_MUTABLE_QUEUE, CGAL_BOOST_PENDING_RELAXED_HEAP };

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

  const value_type& top() const { return mHeap.top(); }

  void pop() { mHeap.pop(); }

  bool empty() const { return mHeap.empty() ; }

  bool contains ( value_type const& v ) { return mHeap.contains(v) ; }

  std::optional<value_type> extract_top()
  {
    std::optional<value_type> r ;
    if ( !empty() )
    {
      value_type v = top();
      pop();
      r = std::optional<value_type>(v) ;
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
struct Modifiable_priority_queue<IndexedType_, Compare_, ID_, CGAL_BOOST_PAIRING_HEAP>
{
  typedef Modifiable_priority_queue Self;

  typedef IndexedType_ IndexedType;
  typedef Compare_     Compare;
  typedef ID_          ID;

  struct Reverse_compare
  {
    Compare c;
    Reverse_compare(const Compare& c):c(c){}
    bool operator() (const IndexedType& a, const IndexedType& b) const
    {
      return !c(a,b);
    }
  };

  // reference time (SMS Iphigenia, keeping 0.05% of edges, all default parameters + Surface_mesh): 12045ms
  // --
  // boost::heap::priority_queue is ummutable and cannot be used
  // --
  // typedef boost::heap::d_ary_heap<IndexedType, boost::heap::arity<2>, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true> > Heap; //(15291ms)
  // typedef boost::heap::d_ary_heap<IndexedType, boost::heap::arity<3>, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true> > Heap; //(14351ms)
  // typedef boost::heap::d_ary_heap<IndexedType, boost::heap::arity<4>, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true> > Heap; //(13869ms)
  // typedef boost::heap::d_ary_heap<IndexedType, boost::heap::arity<5>, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true> > Heap; //(13879ms)
  // typedef boost::heap::d_ary_heap<IndexedType, boost::heap::arity<6>, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true> > Heap; //(13881ms)
  // --
  //typedef boost::heap::binomial_heap<IndexedType, boost::heap::compare<Reverse_compare>> Heap; //(16216ms)
  // --
  // typedef boost::heap::fibonacci_heap<IndexedType, boost::heap::compare<Reverse_compare>> Heap; // (13523ms)
  // --
  typedef boost::heap::pairing_heap<IndexedType, boost::heap::compare<Reverse_compare>> Heap; // (12174ms)
  // --
  // typedef boost::heap::skew_heap<IndexedType, boost::heap::compare<Reverse_compare>, boost::heap::mutable_<true>> Heap; //(17957ms)
  // --

  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;
  typedef typename Heap::handle_type handle_type;

private:
  void reserve_impl(size_type r, std::true_type)
  {
    mHeap.reserve(r);
  }

  void reserve_impl(size_type, std::false_type)
  {}

public:
  Modifiable_priority_queue(size_type largest_ID,
                            const Compare& c = Compare(),
                            const ID& id = ID())
    : mHeap(Reverse_compare(c))
    , mID(id)
    , mHandles(largest_ID, handle_type())
  {
    reserve(largest_ID);
  }

  Modifiable_priority_queue(const Modifiable_priority_queue&) = delete;
  Modifiable_priority_queue& operator=(const Modifiable_priority_queue&) = delete;

public:
  handle_type push(const value_type& v)
  {
    CGAL_precondition(!contains(v));
    auto vid = get(mID, v);
    mHandles[vid] = mHeap.push(v);
    return mHandles[vid];
  }

  handle_type resize_and_push(const value_type& v)
  {
    auto vid = get(mID, v);
    CGAL_precondition(0 <= vid);

    if(vid >= mHandles.size())
    {
      mHandles.resize(vid + 1, handle_type());
//      std::cout << "resize() to " << mHandles.size() << " (capacity = " << mHandles.capacity() << ")" << std::endl;
    }

    CGAL_precondition(!contains(v));
    mHandles[vid] = mHeap.push(v);
    return mHandles[vid];
  }

  void update(const value_type& v)
  {
    CGAL_precondition(contains(v));
    auto vid = get(mID, v);
    mHeap.update(mHandles[vid]);
  }

  void erase(const value_type& v)
  {
    CGAL_precondition(contains(v));
    auto vid = get(mID, v);
    mHeap.erase(mHandles[vid]);
    mHandles[vid] = handle_type();
  }

  const value_type& top() const
  {
    return mHeap.top();
  }

  void pop()
  {
    auto vid = get(mID, top());
    mHeap.pop();
    mHandles[vid] = handle_type();
  }

  size_type size() const { return mHeap.size(); }

  bool empty() const { return mHeap.empty(); }

  void clear()
  {
    std::fill(std::begin(mHandles), std::end(mHandles), handle_type());
    mHeap.clear();
  }

  bool contains(const value_type& v)
  {
    auto vid = get(mID, v);
    CGAL_precondition(0 <= vid && vid < mHandles.size());
    return mHandles[vid] != handle_type();
  }

  bool contains_with_bounds_check(const value_type& v)
  {
    auto vid = get(mID, v);
    if(vid >= mHandles.size())
      return false;

    CGAL_precondition(0 <= vid && vid < mHandles.size());
    return (mHandles[vid] != handle_type());
  }

  std::optional<value_type> extract_top()
  {
    std::optional<value_type> r;
    if(!empty())
    {
      value_type v = top();
      pop();
      r = std::optional<value_type>(v);
    }

    return r;
  }

  value_type top_and_pop()
  {
    CGAL_precondition(!empty());
    value_type v = top();
    pop();
    return v;
  }

  void reserve(size_type r)
  {
    reserve_impl(r, std::integral_constant<bool, Heap::has_reserve>());
  }

private:
  Heap mHeap;
  ID mID;
  std::vector<handle_type> mHandles;
};

} // namespace CGAL

#endif // CGAL_MODIFIABLE_PRIORITY_QUEUE_H
