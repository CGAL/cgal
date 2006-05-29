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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_MODIFIABLE_PRIORITY_QUEUE_H
#define CGAL_MODIFIABLE_PRIORITY_QUEUE_H

#include<CGAL/Multiset.h>

CGAL_BEGIN_NAMESPACE

template <class Type_ 
         ,class Compare_ = CGAL::Compare<Type_>
         ,class Allocator_ = CGAL_ALLOCATOR(int)
          >
class Modifiable_priority_queue
{
  
public:

  typedef Multiset<Type_,Compare_,Allocator_> Heap ;
  
  typedef Modifiable_priority_queue Self;
  
  typedef Compare_ Compare;
  
  typedef typename Heap::value_type      value_type;
  typedef typename Heap::reference       reference;
  typedef typename Heap::const_reference const_reference;
  typedef typename Heap::size_type       size_type;
  typedef typename Heap::iterator        iterator ;
  
public:

  Modifiable_priority_queue() {}
  
  Modifiable_priority_queue( Compare const& c ) : mHeap(c) {}
  
  template<class InputIterator>
  Modifiable_priority_queue( InputIterator first, InputIterator last, Compare const& c )
    :
    mHeap(first,last,c)
  {}  
  
  iterator push ( value_type const& v ) { return mHeap.insert(v) ; }
  
  iterator update ( iterator i ) { value_type v = *i; mHeap.erase(i); return mHeap.insert(v); }
  
  void erase ( iterator i ) { mHeap.erase(i); }

  bool contains ( value_type const& v ) { return mHeap.find(v) != mHeap.end() ; }
    
  value_type top() const { return *mHeap.begin() ; }
  
  void pop() { mHeap.erase(mHeap.begin()); }
  
  bool empty() const { return mHeap.empty() ; }
  
  size_type size() const { return mHeap.size() ; }

  bool is_valid_iterator ( iterator j )
  {
    for( iterator it = mHeap.begin(); it != mHeap.end() ; ++ it )
      if ( j == it )
        return true ;
    return false ;    
  }  

private:

  Heap mHeap ;  
    
} ;

CGAL_END_NAMESPACE

#endif
 
