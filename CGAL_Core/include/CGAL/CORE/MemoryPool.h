/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
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
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/
#ifndef _CORE_MEMORYPOOL_H_
#define _CORE_MEMORYPOOL_H_

#include <CGAL/config.h>
#include <CGAL/tss.h>

#include <new>           // for placement new
#include <cassert>
#include <CGAL/assertions.h>
#include <vector>

namespace CORE {

#define CORE_EXPANSION_SIZE 1024
template< class T, int nObjects = CORE_EXPANSION_SIZE >
class MemoryPool {
private:
   struct Thunk {
      T object;
      Thunk* next;
   };

  typedef MemoryPool<T,nObjects> Self;
public:
   MemoryPool() : head( 0 ) {}

  ~MemoryPool()
  {
    //CGAL_warning_code(
      std::size_t count = 0;
      Thunk* t = head;
      while(t!=0){
        ++count;
        t = t->next;
      }
    //);
    //CGAL_warning_msg(count ==  nObjects * blocks.size(),
    //                 "Cannot delete memory as there are cyclic references");

    if(count ==  nObjects * blocks.size()){
      for(std::size_t i=0; i < blocks.size();i++){
        ::operator delete(blocks[i]);
      }
    }
    // un-commenting the following line can help reproduce on Linux the
    // assertion !blocks.empty() that is sometimes triggered with MSVC
    // or AppleClang
    // blocks.clear();
  }


   void* allocate(std::size_t size);
   void free(void* p);

  // Access the corresponding static global allocator.
  static MemoryPool<T,nObjects>& global_allocator() {
#if defined(CGAL_HAS_THREADS) // use the C++11 implementation
    static thread_local Self memPool;
#else // not CGAL_HAS_THREADS
    static Self memPool;
#endif // not CGAL_HAS_THREADS
    return memPool;
  }

private:
   Thunk* head; // next available block in the pool
  std::vector<void*> blocks;
};

template< class T, int nObjects >
void* MemoryPool< T, nObjects >::allocate(std::size_t) {
   if ( head == 0 ) { // if no more memory in the pool
      const int last = nObjects - 1;

      // use the global operator new to allocate a block for the pool
      Thunk* pool = reinterpret_cast<Thunk*>(
         ::operator new(nObjects * sizeof(Thunk)));

      blocks.push_back(pool);
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
   if(blocks.empty()){
     std::cerr << typeid(T).name() << std::endl;
   }
   CGAL_assertion (! blocks.empty());

   // recycle the object memory, by putting it back into the chain
   reinterpret_cast<Thunk*>(t)->next = head;
   head = reinterpret_cast<Thunk*>(t);
}

} //namespace CORE

#endif // _CORE_MEMORYPOOL_H_
