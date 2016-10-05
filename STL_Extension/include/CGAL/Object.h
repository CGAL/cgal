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

#include <CGAL/config.h>
#include <CGAL/assertions.h>

#include <typeinfo>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

namespace CGAL {

class Object
{
    boost::shared_ptr<boost::any> obj;
  
    // returns an any pointer from a variant
    struct Any_from_variant : public boost::static_visitor<boost::any*> {
      template<typename T>
      boost::any* operator()(const T& t) const { 
        return new boost::any(t);
      }  
    };

    template<class T>
    friend const T* object_cast(const Object * o);

    template<class T>
    friend T object_cast(const Object & o);

    typedef void (Object::*bool_type)() const;
    void this_type_does_not_support_comparisons() const {}
  public:

    struct private_tag{};

    Object() : obj() { }

#ifndef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
    template <class T>
    Object(T && t, private_tag) : obj(new boost::any(std::forward<T>(t))) { }
#else
    template <class T>
    Object(const T&t, private_tag) : obj(new boost::any(t)) { }
#endif

    // implicit constructor from optionals containing variants
    template<BOOST_VARIANT_ENUM_PARAMS(typename T)>
    Object(const boost::optional< boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) > >& t) 
      : obj( t ? boost::apply_visitor(Any_from_variant(), *t) : NULL) { }
  
    // implicit constructor from  variants
    template<BOOST_VARIANT_ENUM_PARAMS(typename T)>
    Object(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& v) 
      : obj(boost::apply_visitor(Any_from_variant(), v)) { }

    template <class T>
    bool assign(T &t) const
    {
      if(obj) {
        #ifdef CGAL_USE_ANY_BAD_CAST
        try {
          t = boost::any_cast<T>(*obj);
          return true;
        } catch(...) {
          return false;
        }
        #else
        const T* res =  boost::any_cast<T>(&(*obj));
        if (res){
          t=*res;
          return true;
        }
        return false;
        #endif
      } else {
        return false;
      }
    }

    bool
    empty() const
    {
      return !obj;
    }

    // is_empty() is kept for backward compatibility.
    // empty() was introduced for consistency with e.g. std::vector::empty().
    bool
    is_empty() const
    {
	return empty();
    }

    // safe-bool conversion
    operator bool_type() const {
      return empty() == false ? &Object::this_type_does_not_support_comparisons : 0;
    }


    template <class T>
    bool is() const
    {
      return obj && boost::any_cast<T>(obj.get());
    }

    const std::type_info & type() const
    {
      if(obj)
        return obj->type();
      else
        return typeid(void);
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
  if(o->obj)
    return boost::any_cast<T>((o->obj).get());
  else
    return NULL;
}

template <class T>
inline
T object_cast(const Object & o)
{
  if(!o.obj)
    throw Bad_object_cast();
  
  const T * result = boost::any_cast<T>((o.obj).get());
  if (!result)
    throw Bad_object_cast();
  return *result;
}

} //namespace CGAL

#endif // CGAL_OBJECT_H
