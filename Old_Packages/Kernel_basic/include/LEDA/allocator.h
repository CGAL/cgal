// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, December 10
//
// file          : ?/include/LEDA/allocator.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// author(s)     : ?
// coordinator   : ?
//
// ======================================================================
#ifndef LEDA_ALLOCATOR_H
#define LEDA_ALLOCATOR_H

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 390009
#include <LEDA/REDEFINE_NAMES.h>
#endif


#include <LEDA/basic.h>

// the following piece of code is programmed according to
// the C++ standard clause 20.4.1

template <class T> class leda_allocator;

// specialize for void:
template <> class leda_allocator<void> {
public:
  typedef void*       pointer;
  typedef const void* const_pointer;
  //  reference-to-void members are impossible.
  typedef void        value_type;
  template <class U> struct rebind { typedef leda_allocator<U> other; };
};

/*{\Manpage {leda_allocator} {T} {Memory Allocator} {A}}*/

template <class T>
class leda_allocator {

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is
a memory allocator according to the \CC standard. |\Mname| is the
standard compliant interface to the LEDA memory management.}*/

public:
/*{\Mtypes 5}*/
/*{\Mtext Local types are |size_type|, |difference_type|, |value_type|,
|pointer|, |reference|, |const_pointer|, and |const_reference|.}*/

  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;
  typedef T           value_type;
  typedef T*          pointer;
  typedef const T*    const_pointer;
  typedef T&          reference;
  typedef const T&    const_reference;

  template <class T1> class rebind { public:
  /*{\Mtypemember allows the construction of a derived allocator:\\
     |\Mname::template rebind<T1>::other|\\ is the type
     |leda_allocator<T1>|. }*/
    typedef leda_allocator<T1> other;
  };


/*{\Mcreation 4.5}*/

  leda_allocator() {}
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname|. }*/

  template <class TO> 
  leda_allocator(const leda_allocator<TO>&) {}
  ~leda_allocator() {}

/*{\Moperations 3 3 }*/

pointer allocate(size_type n, const_pointer = 0)
/*{\Mop returns a pointer to a newly allocated memory range of size
        |n * sizeof(T)|.}*/
{ return 0 == n ? 0 : 
         (T*) std_memory_mgr.allocate_bytes( n * sizeof(T) ); }

void deallocate(pointer p, size_type n)
/*{\Mop deallocates a memory range of |n * sizeof(T)| starting
        at |p|. \precond the memory range was obtained via |allocate(n)|.}*/
{ std_memory_mgr.deallocate_bytes(p , n * sizeof(T)); }

pointer address(reference r)
/*{\Mop returns |&r|.}*/
{ return &r; }

const_pointer address(const_reference r)
/*{\Mop returns |&r|.}*/
{ return &r; }

void construct(pointer p, const_reference r)
/*{\Mop makes an inplace new |new( (void*)p ) T(r)|.}*/
{ new(p) value_type(r); }

void destroy(pointer p)
/*{\Mop destroys the object referenced via |p| by calling |p->~T()|.}*/
{ p->~T(); }

size_type max_size() const { return std_memory_mgr.max_size(); }
/*{\Mop the largest value |n| for which the call |allocate(n,0)| 
    might succeed.}*/

};

/*{\Mimplementation Note that the above class template uses all kinds
of modern compiler technology like member templates, partial specialization
etc. It runs only on a subset of LEDA's general supported platforms like
|g++ > 2.95|, |SGI CC > 7.3|.}*/

#if LEDA_ROOT_INCL_ID == 390009
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif // LEDA_ALLOCATOR_H

