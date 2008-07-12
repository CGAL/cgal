// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//                 Sylvain Pion
// ======================================================================

#ifndef CGAL_OBJECT_HANDLE_H
#define CGAL_OBJECT_HANDLE_H

#include <CGAL/Handle_for_virtual.h>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

template <class T>
class Handle_wrapper : public Ref_counted_virtual
{
  public:
    Handle_wrapper(const T& object) : _object(object) {}
    Handle_wrapper() {}
    operator T() const { return _object; }
    ~Handle_wrapper() {}
  private:
  T _object;
};


class Object_handle
  : public Handle_for_virtual<Ref_counted_virtual>
{
  struct empty{};
  typedef Handle_for_virtual<Ref_counted_virtual> base;
public:
  Object_handle()
  { initialize_with(Handle_wrapper<empty>()); }

  template <class T>
  Object_handle(const T&t)
  { initialize_with(Handle_wrapper<T>(t)); }

  template <class T>
  bool assign(T &t) const
  {

      const Handle_wrapper<T> *wp = 
      	dynamic_cast<const Handle_wrapper<T> *>(Ptr());
      if ( wp == static_cast<Handle_wrapper<T> *>(0) )
      	return false;
      t = *(wp);

    return true;
  }

  bool is_empty() const
  { empty E; return assign(E); }
  bool operator==(Nullptr_t n) const
  { CGAL_assertion(n == 0); return is_empty(); }
  bool operator!=(Nullptr_t n) const
  { CGAL_assertion(n == 0); return !is_empty(); }

};

template <class T>
inline bool assign(T& t, const Object_handle& o)
{ return o.assign(t); }

CGAL_END_NAMESPACE

#endif // CGAL_OBJECT_HANDLE_H
