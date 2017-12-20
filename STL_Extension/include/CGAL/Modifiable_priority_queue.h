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
#include <CGAL/internal/boost/mutable_queue.hpp>


namespace CGAL {
  namespace internal {
template <class IndexedType, 
          class RandomAccessContainer = std::vector<IndexedType>, 
          class Comp = std::less<typename RandomAccessContainer::value_type>,
          class ID = ::boost::identity_property_map >
class mutable_queue_with_remove : public internal::boost_::mutable_queue<IndexedType,RandomAccessContainer,Comp,ID>
{
  typedef internal::boost_::mutable_queue<IndexedType,RandomAccessContainer,Comp,ID> Base;
public:
  typedef typename Base::size_type size_type;
  typedef typename Base::Node Node;

  mutable_queue_with_remove(size_type n, const Comp& x=Comp(), const ID& _id=ID()) : Base(n,x,_id,true) 
  {}

  void remove(const IndexedType& x){
    //first place element at the top
    size_type current_pos = this->index_array[ get(this->id, x) ];
    this->c[current_pos] = x;

    Node node(this->c.begin(), this->c.end(), this->c.begin()+current_pos, this->id);
    while (node.has_parent())
      node.swap(node.parent(), this->index_array);
    //then pop it
    this->pop();
  }
  
  bool contains(const IndexedType& x) const {
    return this->index_array[ get(this->id, x) ] !=this->index_array.size();
  }
};

} } //namespace CGAL::internal
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
  typedef  internal::mutable_queue_with_remove<IndexedType,std::vector<IndexedType>,Compare,ID> Heap;
  #endif //CGAL_SURFACE_MESH_SIMPLIFICATION_USE_RELAXED_HEAP
  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;
  
  typedef bool handle ;
  
public:

  Modifiable_priority_queue( size_type largest_ID, Compare const& c, ID const& id ) : mHeap(largest_ID,c,id) {}
  
  handle push ( value_type const& v ) { mHeap.push(v) ; return handle(true) ; }
  
  handle update ( value_type const& v, handle h ) { mHeap.update(v); return h ; }
  
  handle erase ( value_type const& v, handle  ) { mHeap.remove(v); return null_handle() ; }
  handle erase ( value_type const& v  ) { mHeap.remove(v); return null_handle() ; }

  value_type top() const { return mHeap.top() ; }
  
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
  
  static handle null_handle() { return handle(false); }
  
private:

  Heap mHeap ;  
    
} ;

} //namespace CGAL

#endif
 
