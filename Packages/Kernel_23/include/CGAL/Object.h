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

#include <CGAL/basic.h>
#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

template <class T>
class Wrapper : public Ref_counted_virtual
{
  public:
    Wrapper(const T& object) : _object(object) {}

    Wrapper() {}

    operator T() const { return _object; }

    ~Wrapper() {}

  private:
    T         _object;
};

class Object
  : public Handle_for_virtual<Ref_counted_virtual>
{
    struct empty{};
    typedef Handle_for_virtual<Ref_counted_virtual> base;

  public:

    struct private_tag{};

    Object()
    {
	initialize_with(Wrapper<empty>());
    }

    template <class T>
    Object(const T&t, private_tag)
    {
	initialize_with(Wrapper<T>(t));
    }

    template <class T>
    bool assign(T &t) const
    {
#ifdef _MSC_VER
      try {
#endif
        const Wrapper<T> *wp = dynamic_cast<const Wrapper<T> *>(Ptr());
        if (wp == NULL)
            return false;
        t = *wp;
#ifdef _MSC_VER
      }
      catch (...) {
          std::cerr << "ERROR : YOUR COMPILER MUST SUPPORT RTTI" << std::endl;
          abort();
      }
#endif
      return true;
    }

    bool
    is_empty() const
    {
	empty E;
	return assign(E);
    }
};


template <class T>
inline
Object
make_object(const T& t)
{
    return Object(t, Object::private_tag());
}

template <class T>
inline
bool
assign(T& t, const Object& o)
{
    return o.assign(t);
}

CGAL_END_NAMESPACE

#endif // CGAL_OBJECT_H
