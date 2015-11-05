/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: MemoryPool.h
 * Synopsis:
 *      a memory pool template class.
 * 
 * Written by 
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/
#ifndef _CORE_MEMORYPOOL_H_
#define _CORE_MEMORYPOOL_H_

#include <new>           // for placement new
#include <cassert>
#include <CGAL/assertions.h>

namespace CORE { 

#define CORE_EXPANSION_SIZE 1024
template< class T, int nObjects = CORE_EXPANSION_SIZE >
class MemoryPool {
public:
   MemoryPool() : head( 0 ) {}

   void* allocate(std::size_t size);
   void free(void* p);

  // Access the corresponding static global allocator.
  static MemoryPool<T>& global_allocator() {
    return memPool;
  }
  
private:
   struct Thunk { 
      T object;
      Thunk* next;
   };

private:
   Thunk* head; // next available block in the pool

private:
  // Static global allocator.
  static MemoryPool<T, nObjects> memPool;   
};

template <class T, int nObjects >
MemoryPool<T, nObjects> MemoryPool<T, nObjects>::memPool;

template< class T, int nObjects >
void* MemoryPool< T, nObjects >::allocate(std::size_t) {
   if ( head == 0 ) { // if no more memory in the pool
      const int last = nObjects - 1;

      // use the global operator new to allocate a block for the pool
      Thunk* pool = reinterpret_cast<Thunk*>(
	 ::operator new(nObjects * sizeof(Thunk)));

      // initialize the chain (one-directional linked list)
      head = pool;
      for (int i = 0; i < last; ++i ) {
	 pool[i].next = &pool[i+1];
      }
      pool[last].next = 0;
   }

   // set head to point to the next object in the list
   Thunk* currentThunk = head;
   head = currentThunk->next;

   return currentThunk;
}

template< class T, int nObjects >
void MemoryPool< T, nObjects >::free(void* t) {
   CGAL_assertion(t != 0);     
   if (t == 0) return; // for safety

   // recycle the object memory, by putting it back into the chain
   reinterpret_cast<Thunk*>(t)->next = head;
   head = reinterpret_cast<Thunk*>(t);
}

} //namespace CORE

#endif // _CORE_MEMORYPOOL_H_
