// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// source        : No_op_allocator.fw
// file          : No_op_allocator.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.17
// revision_date : 10 Dec 2000 
// author(s)     : Lutz Kettner
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_NO_OP_ALLOCATOR_H
#define CGAL_NO_OP_ALLOCATOR_H

template <class T>
class No_op_allocator
{
 public:
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;
  typedef T           value_type;
  typedef T*          pointer;
  typedef const T*    const_pointer;
  typedef T&          reference;
  typedef const T&    const_reference;


  No_op_allocator() : Don_t_call(true) {}

  pointer
  allocate(size_type n, const_pointer = 0) { return static_cast<pointer>(0); }

  void
  deallocate(pointer p, size_type n) {}

  pointer
  address(reference r)
  { CGAL_precondition ( !Don_t_call ); return static_cast<pointer>(0); }

  const_pointer
  address(const_reference r)
  { CGAL_precondition ( !Don_t_call ); return static_cast<const_pointer>(0); }

  void
  construct(pointer p, const_reference r)
  { CGAL_precondition ( !Don_t_call ); }
  // It is the user responsibility to construct
  // the element, e.g. using new

  void
  destroy(pointer p)
  { delete p; }
  // p should have an appropriate type
  // That's why we need the parameterization
  // T should have a virtual destructor

  size_type
  max_size() const
  { return 0; }

 private:
  bool  Don_t_call;
};

#endif // CGAL_NO_OP_ALLOCATOR_H
