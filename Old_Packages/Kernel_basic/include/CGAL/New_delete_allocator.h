// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : New_delete_allocator.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_NEW_DELETE_ALLOCATOR_H
#define CGAL_NEW_DELETE_ALLOCATOR_H

template <class T>
class New_delete_allocator
{
 public:
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;
  typedef T           value_type;
  typedef T*          pointer;
  typedef const T*    const_pointer;
  typedef T&          reference;
  typedef const T&    const_reference;

  template <class U> struct rebind { typedef New_delete_allocator<U> other; };

  New_delete_allocator() {}

  pointer
  allocate(size_type, const_pointer = 0)
  { return static_cast<pointer>(0); }

  void
  deallocate(pointer, size_type) {}

  pointer
  address(reference)
  { return static_cast<pointer>(0); }

  const_pointer
  address(const_reference)
  { return static_cast<const_pointer>(0); }

  template <class U>
  void
  construct(U*& ptr, const U& ref)
  { ptr = new U(ref); }

  void
  destroy(pointer p)
  { delete p; }
  // p should have an appropriate type
  // That's why we need the parameterization
  // T should have a virtual destructor

  size_type
  max_size() const
  { return 0; }
};

#endif // CGAL_NEW_DELETE_ALLOCATOR_H
