// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
