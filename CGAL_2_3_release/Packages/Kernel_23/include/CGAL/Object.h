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
// release_date  : 
//
// file          : include/CGAL/Object.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan.Schirra  <Stefan.Schirra@mpi-sb.mpg.de>
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//
// coordinator   : MPI Saarbruecken, Germany 
//
// ======================================================================

#ifndef CGAL_OBJECT_H
#define CGAL_OBJECT_H

#ifdef CGAL_CFG_NO_DYNAMIC_CAST
#error fatal error: dynamic cast not supported
#endif

#include <CGAL/Handle_for.h>

namespace CGAL {

class Object;
class Object_base;
template <class T> class Wrapper;


class Object_base : public Ref_counted  
{
  public:
    virtual   ~Object_base() {}
};


template <class T>
class Wrapper : public Object_base 
{
  public:
    Wrapper(const T& object) : _object(object) {}

    Wrapper() {}

    operator T() { return _object; }

    virtual   ~Wrapper() {}

  private:
    T         _object;
};


class Object
{
  public:
    Object() : ptr( static_cast<Object_base*>(0) ) {}

    Object(Object_base *base) 
    { 
      ptr = base; 
      CGAL_kernel_assertion( !ptr || (ptr->count == 1));
    }

    Object(const Object& o) : ptr(o.ptr)
    { if (ptr) ptr->count++; }

    ~Object()
    { if (ptr && (--ptr->count == 0)) { delete ptr; } }

    Object&       
    operator=(const Object& o)
    {
      if (o.ptr) o.ptr->count++;
      if (ptr && (--ptr->count == 0)) { delete ptr; }
      ptr = o.ptr;
      return *this;
    }

    bool          
    is_empty() const { return ptr == static_cast<Object_base*>(0); }

    template <class T>
    friend bool assign(T& t, const Object& o);

  protected:
    Object_base*  ptr;
};


template <class T>
Object
make_object(const T& t)
{ return Object(new Wrapper< T >(t)); }


template <class T>
bool
assign(T& t, const Object& o)
{
# ifdef CGAL_CFG_DYNAMIC_CAST_BUG
  Wrapper<T>   instantiate_it;
# endif // CGAL_CFG_DYNAMIC_CAST_BUG
  Wrapper<T>*  wp = dynamic_cast<Wrapper<T>*>(o.ptr);
  if ( wp == static_cast<Wrapper<T>*>(0) ) { return false; }
  t = *(wp);
  return true;
}

} // namespace CGAL

#endif // CGAL_OBJECT_H

