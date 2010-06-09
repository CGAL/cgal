// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_MODIFIABLE_PRIORITY_QUEUE_H
#define CGAL_MODIFIABLE_PRIORITY_QUEUE_H

#include <climits> // Neeeded by the following Boost header for CHAR_BIT.
#include <boost/pending/relaxed_heap.hpp>

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
  
  typedef boost::relaxed_heap<IndexedType,Compare,ID> Heap ;
  
  typedef typename Heap::value_type value_type;
  typedef typename Heap::size_type  size_type;
  
  typedef bool handle ;
  
public:

  Modifiable_priority_queue( size_type largest_ID, Compare const& c, ID const& id ) : mHeap(largest_ID,c,id) {}
  
  handle push ( value_type const& v ) { mHeap.push(v) ; return handle(true) ; }
  
  handle update ( value_type const& v, handle h ) { mHeap.update(v); return h ; }
  
  handle erase ( value_type const& v, handle  ) { mHeap.remove(v); return null_handle() ; }

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
 
