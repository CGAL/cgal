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
// author(s)     : Stefan Schirra
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//                 Sylvain Pion
//
// coordinator   : MPI Saarbruecken, Germany 
//
// ======================================================================

#ifndef CGAL_OBJECT_H
#define CGAL_OBJECT_H

#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

class Object_base : public Ref_counted_virtual
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

    operator T() const { return _object; }

    virtual   ~Wrapper() {}

  private:
    T         _object;
};


class Object
  : public Handle_for_virtual<Object_base>
{
    struct empty{};
    typedef Handle_for_virtual<Object_base> base;

  public:

    template <class T>
    friend bool assign(T& t, const Object& o);

    Object()
    {
	initialize_with(Wrapper<empty>());
    }

    Object(const Object &o)
	: base(o) {}

    template <class T>
    Object(const T&t)
    {
	initialize_with(Wrapper<T>(t));
    }

    bool
    is_empty() const
    {
	empty E;
	return assign(E, *this);
    }
};


template <class T>
Object
make_object(const T& t)
{
    return Object(t);
}


template <class T>
bool
assign(T& t, const Object& o)
{
#ifdef _MSC_VER
  try {
#endif
    const Wrapper<T> *wp = dynamic_cast<const Wrapper<T> *>(o.Ptr());
    if ( wp == static_cast<Wrapper<T> *>(0) )
        return false;
    t = *(wp);
#ifdef _MSC_VER
  }
  catch (...) {
      std::cerr << "ERROR : YOUR COMPILER MUST SUPPORT RTTI" << std::endl;
      abort();
  }
#endif
  return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_OBJECT_H
