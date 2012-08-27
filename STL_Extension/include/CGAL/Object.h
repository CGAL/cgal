// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
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

namespace CGAL {

template <class T>
class Wrapper : public Ref_counted_virtual
{
    Wrapper(const Wrapper&); // deleted to make sure we don't make useless copies
  public:

    Wrapper(const T& object) : _object(object) {}

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    Wrapper(T && object) : _object(std::move(object)) {}
#endif

    ~Wrapper() {}

    const T& get() const { return _object; }

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
    struct Empty {};
    typedef Handle_for_virtual<Ref_counted_virtual> base;

  public:

    struct private_tag{};

    Object()
    {
	typedef Wrapper<Empty>  Wrap;
        ptr = new Wrap(Empty());
    }

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    template <class T>
    Object(T && t, private_tag)
    {
	typedef Wrapper< typename std::remove_cv< typename std::remove_reference<T>::type >::type >  Wrap;
        ptr = new Wrap(std::forward<T>(t));
    }
#else
    template <class T>
    Object(const T&t, private_tag)
    {
	typedef Wrapper<T>  Wrap;
        ptr = new Wrap(t);
    }
#endif

    template <class T>
    bool assign(T &t) const
    {
#ifdef _MSC_VER
      try {
#endif
        const Wrapper<T> *wp = dynamic_cast<const Wrapper<T> *>(Ptr());
        if (wp == NULL)
            return false;
        t = wp->get();
#ifdef _MSC_VER
      }
      catch (...) {
          CGAL_error_msg("Your compiler must support Run-Time Type Information (RTTI)");
      }
#endif
      return true;
    }

    bool
    empty() const
    {
	Empty E;
	return assign(E);
    }

    // is_empty() is kept for backward compatibility.
    // empty() was introduced for consistency with e.g. std::vector::empty().
    bool
    is_empty() const
    {
	return empty();
    }

    template <class T>
    bool is() const
    {
        return NULL != dynamic_cast<const Wrapper<T> *>(Ptr());
    }

    const std::type_info & type() const
    {
        return empty() ? typeid(void) : Ptr()->type();
    }

#ifndef CGAL_NO_DEPRECATED_CODE
    // The comparisons with NULL are only there for Nef...
    bool operator==(Nullptr_t CGAL_assertion_code(n)) const
    { CGAL_assertion(n == 0); return empty(); }
    bool operator!=(Nullptr_t CGAL_assertion_code(n)) const
    { CGAL_assertion(n == 0); return !empty(); }
#endif // CGAL_NO_DEPRECATED_CODE

};


#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
template <class T>
inline
Object
make_object(T && t)
{
    return Object(std::forward<T>(t), Object::private_tag());
}
#else
template <class T>
inline
Object
make_object(const T& t)
{
    return Object(t, Object::private_tag());
}
#endif

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


template <class T>
inline
const T * object_cast(const Object * o)
{
    const Wrapper<T> *wp = dynamic_cast<const Wrapper<T> *>(o->Ptr());
    if (wp == NULL)
        return NULL;
    return static_cast<const T*>(wp->object_ptr());
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

} //namespace CGAL

#endif // CGAL_OBJECT_H
