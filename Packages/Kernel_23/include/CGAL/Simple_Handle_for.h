// Copyright (c) 2001,2002,2003  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_SIMPLE_HANDLE_FOR_H
#define CGAL_SIMPLE_HANDLE_FOR_H

CGAL_BEGIN_NAMESPACE

template < class Stored >
class Simple_Handle_for
{
public:

    typedef Stored element_type;

    Simple_Handle_for()
    {}

    Simple_Handle_for(const Stored& rc)
	: _s(rc) {}

/* I comment this one for now, since it's preventing the automatic conversions
   to take place.  We'll see if it's a problem later.
    template < typename T1 >
    Simple_handle_for(const T1& t1)
      : _s(t1) {}
*/

    template < typename T1, typename T2 >
    Simple_Handle_for(const T1& t1, const T2& t2)
      : _s(t1, t2) {}

    template < typename T1, typename T2, typename T3 >
    Simple_Handle_for(const T1& t1, const T2& t2, const T3& t3)
      : _s(t1, t2, t3) {}

    template < typename T1, typename T2, typename T3, typename T4 >
    Simple_Handle_for(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
      : _s(t1, t2, t3, t4) {}

    Simple_Handle_for&
    operator=(const Stored& rc)
    {
      _s = rc;
      return *this;
    }

    void
    initialize_with(const Stored& rc)
    {
      _s = rc;
    }

    long int
    id() const
    { return reinterpret_cast<long int>(&_s); }

    bool
    identical(const Simple_Handle_for& h) const
    { return id() == h.id(); } // Or should it always return false ?

    const Stored * Ptr() const
    { return &_s; }

    Stored * Ptr()
    { return &_s; }

    const Stored * ptr() const
    { return &_s; }

    Stored * ptr()
    { return &_s; }

    bool
    is_shared() const
    {
	return false;
    }

    void
    swap(Simple_Handle_for &h)
    {
        using std::swap;
        swap(_s, h._s);
    }

private:
    Stored _s;
};

template < class Stored >
inline
void
swap(Simple_Handle_for<Stored> &h1, Simple_Handle_for<Stored> &h2)
{
    h1.swap(h2);
}

template < class Stored >
inline
bool
identical(const Simple_Handle_for<Stored> &h1,
          const Simple_Handle_for<Stored> &h2)
{
    return h1.identical(h2);
}

template <class T>
inline
const T&
get(const Simple_Handle_for<T> &h)
{
    return *(h.Ptr());
}

CGAL_END_NAMESPACE

#endif // CGAL_SIMPLE_HANDLE_FOR_H
