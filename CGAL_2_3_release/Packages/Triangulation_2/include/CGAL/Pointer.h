// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Pointer.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers, Mariette Yvinec, Sylvain Pion
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_POINTER_H
#define CGAL_POINTER_H

#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE 

template <class T>
struct Pointer
{
  typedef T value_type;

  //Pointer(const T* p = NULL) : _pointer((T*)p) {}
  Pointer( T* p = NULL) : _pointer((T*)p) {}

  //Pointer& operator=(const T* p)
  Pointer& operator=( T* p)
    {
      ptr() = p;
      return *this;
    }

  Pointer& operator=(const Pointer& p)
    {
      ptr() = p.ptr();
      return *this;
    }

  T& operator*() const { return *ptr(); }
  T* operator->() const { return ptr(); }

  void clear() { ptr() = NULL; }

  void Delete()
    {
      delete ptr();
      clear();
    }

  bool is_null() const { return ptr() == NULL; }

  bool operator==(const Pointer& p) const { return ptr() == p.ptr(); }
  bool operator!=(const Pointer& p) const { return !(*this == p); }

  bool operator==(CGAL_NULL_TYPE CGAL_assertion_code(n) ) const
    {
      CGAL_assertion(n == 0);
      return ptr() == NULL;
    }

  bool operator!=(CGAL_NULL_TYPE n) const { return !(*this == n); }
 
  T*& ptr()       { return _pointer; }
  T*  ptr() const { return _pointer; }

  bool operator<(const Pointer& p) const { return ptr() < p.ptr();}

private:
  T* _pointer;
};

CGAL_END_NAMESPACE

#endif // CGAL_POINTER_H
