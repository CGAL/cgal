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
// release_date  : 
// 
// file          : No_op_allocator.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
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
  allocate(size_type, const_pointer = 0)
  { return static_cast<pointer>(0); }

  void
  deallocate(pointer, size_type) {}

  pointer
  address(reference)
  { CGAL_precondition ( !Don_t_call ); return static_cast<pointer>(0); }

  const_pointer
  address(const_reference)
  { CGAL_precondition ( !Don_t_call ); return static_cast<const_pointer>(0); }

  void
  construct(pointer, const_reference)
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
