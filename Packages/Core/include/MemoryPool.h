/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: MemoryPool.h
 * Synopsis:
 *      a memory pool template class.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_MEMORYPOOL_H
#define CORE_MEMORYPOOL_H

#ifdef CORE_DISABLE_MEMORY_POOL
#define CORE_MEMORY(T)
#else

#include <cstdlib>
#include <typeinfo>

CORE_BEGIN_NAMESPACE

#define CORE_EXPANSION_SIZE 1024

template<class T>
class MemoryPool {
public:
  MemoryPool();
  ~MemoryPool();
  
  // allocate a T element from the free list.
  void* allocate(size_t) {
    if (next == NULL)
      expandFreeList();

    MemoryPool<T> *head = next;
    next = head->next;
    
    nCount --;

    return head;
  }

  // return a T element to the free list.
  void free(void* doomed) {
    if (doomed == NULL)
      return;

    MemoryPool<T> *head = static_cast<MemoryPool<T> *>(doomed);
   
    head->next = next;
    next = head;

    nCount ++;
    //if (nCount % 2048 == 0) {
    //  std::cout << nCount << " of " << typeid(T).name() 
    //       << "(size=" << sizeof(T) << ")" << ", total=" 
    //       << (nCount>>10)*sizeof(T) << "KB" << std::endl;
    //  releaseFreeList(1024);
    //}
  }

  // Access the corresponding static global allocator.
  static MemoryPool<T>& global_allocator()
  {
    return memPool;
  }

private:
  // next element on the free list.
  MemoryPool<T> *next;
  
  // expand free list.
  void expandFreeList(int howMany = CORE_EXPANSION_SIZE);

  // release free list.
  void releaseFreeList(int howMany);

  // number of the free element.
  int nCount;

  // Static global allocator.
  static MemoryPool<T> memPool;
};

template <class T>
MemoryPool<T> MemoryPool<T>::memPool;

template<class T>
MemoryPool<T>::MemoryPool() : next(NULL), nCount(0) {}

template<class T>
MemoryPool<T>::~MemoryPool() {
  releaseFreeList(nCount);
}

template<class T>
void MemoryPool<T>::releaseFreeList(int howMany) {
  int i;
  MemoryPool<T> *nextPtr = next;  
  for (i=0; (nextPtr != NULL && i < howMany); i++) {
    nextPtr = next;
    next = nextPtr->next;
    ::delete[] reinterpret_cast<char *>(nextPtr);
    //::delete[] reinterpret_cast<T *>(nextPtr);
    //::delete nextPtr;
  }
  nCount -= i; /* in case failure */
}

template<class T>
void MemoryPool<T>::expandFreeList(int howMany) {

  size_t size = sizeof(T);
  if ( size < sizeof(MemoryPool<T> *) )
    size = sizeof(MemoryPool<T> *);

  char* p = ::new char[size];
  if (p == NULL)
    std::cerr << "Out of Memory!!!" << std::endl;

  MemoryPool<T> *runner = reinterpret_cast<MemoryPool<T> *>(p);
  next = runner;
  
  for (int i=0; i<howMany-1; i++) {
    p = ::new char[size];
    if (p == NULL)
      std::cerr << "Out of Memory!!!" << std::endl;

    runner->next = reinterpret_cast<MemoryPool<T> *>(p);
    runner = runner->next;
  }

  runner->next = NULL;
  
  nCount += howMany;

  /* Below implementation will be faster, but if we use it, there is no way to
     free allocated memory except when the program terminated.
  */
  /*
  char* p = new char[size*howMany];
  if (p == NULL)
    std::cerr << "Out of Memory!!!" << std::endl;

  MemoryPool<T> *runner = reinterpret_cast<MemoryPool<T> *>(p);
  next = runner;

  for (int i=0; i<howMany-1; i++) {
    p += size;
    runner->next = reinterpret_cast<MemoryPool<T> *>(p);
    runner = runner->next;
  }

  runner->next = NULL;
  */

}

CORE_END_NAMESPACE

// You can put the following macro in (the public part of) the body of a class
// to make it use the memory pool.

#define CORE_MEMORY(T)                                                 \
  void *operator new( size_t size)                                     \
  { return MemoryPool<T>::global_allocator().allocate(size); }         \
  void operator delete( void *p, size_t )                              \
  { MemoryPool<T>::global_allocator().free(p); }

#endif
#endif
