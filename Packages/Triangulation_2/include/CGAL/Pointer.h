// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Olivier Devillers, Mariette Yvinec, Sylvain Pion

#ifndef CGAL_POINTER_H
#define CGAL_POINTER_H

#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE 

template <class T>
struct Pointer
{
  typedef void     iterator_category;
  typedef T        value_type;
  typedef void     difference_type;
  typedef T*       pointer;
  typedef T&       reference;
  

  //Pointer(const T* p = NULL) : _pointer((T*)p) {}
  //Pointer( T* p = NULL) : _pointer((T*)p) {}
  Pointer() : _pointer(NULL) {}
  Pointer( const T* p) : _pointer(const_cast<T*>(p)) {}

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
  bool operator>(const Pointer& p) const { return ptr() > p.ptr();}

private:
  T* _pointer;
};

CGAL_END_NAMESPACE

#endif // CGAL_POINTER_H
