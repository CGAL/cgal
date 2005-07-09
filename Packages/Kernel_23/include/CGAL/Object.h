// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Stefan Schirra
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//                 Sylvain Pion

#ifndef CGAL_OBJECT_H
#define CGAL_OBJECT_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for_virtual.h>

#include <typeinfo>

CGAL_BEGIN_NAMESPACE

template <class T>
class Wrapper : public Ref_counted_virtual
{
  public:
    Wrapper(const T& object) : _object(object) {}

    Wrapper() {}

    operator T() const { return _object; }

    ~Wrapper() {}

    virtual const std::type_info & type() const
    {
        return typeid(T);
    }

    virtual const void * object_ptr() const
    {
        return & _object;
    }

  private:
    T         _object;
};

class Object
  : public Handle_for_virtual<Ref_counted_virtual>
{
    struct empty {};
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

    const std::type_info & type() const
    {
        return is_empty() ? typeid(void) : Ptr()->type();
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


struct Bad_object_cast
  : public std::bad_cast
{
    virtual const char * what() const throw()
    {
        return "CGAL::bad_object_cast: "
               "failed conversion using CGAL::object_cast";
    }
};

/*
template <class T>
inline
T * object_cast(Object * o)
{
    return o && o->type() == typeid(T) ? o->object_ptr() : NULL;
}
*/

template <class T>
inline
const T * object_cast(const Object * o)
{
    const Wrapper<T> *wp = dynamic_cast<const Wrapper<T> *>(o->Ptr());
    if (wp == NULL)
        return NULL;
    return static_cast<const T*>(wp->object_ptr());
//    return o && o->type() == typeid(T)
//         ? static_cast<const T*>(o->object_ptr())
//         : NULL;
}

template <class T>
inline
T object_cast(const Object & o)
{
    const T * result = object_cast<T>(&o);
    if (!result)
        throw Bad_object_cast();
    return *result;
}

CGAL_END_NAMESPACE

#endif // CGAL_OBJECT_H
