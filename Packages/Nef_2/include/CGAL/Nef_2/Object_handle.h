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
// release       : $CGAL_Revision: CGAL-2.4-I-33 $
// release_date  : $CGAL_Date: 2001/12/04 $
//
// file          : include/CGAL/Object_handle.h
// package       : Nef_2
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>    
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>     
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//                 Sylvain Pion
//
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
  bool operator==(CGAL_NULL_TYPE n) const
  { CGAL_assertion(n == 0); return is_empty(); }
  bool operator!=(CGAL_NULL_TYPE n) const
  { CGAL_assertion(n == 0); return !is_empty(); }

};

template <class T>
inline bool assign(T& t, const Object_handle& o)
{ return o.assign(t); }

CGAL_END_NAMESPACE

#endif // CGAL_OBJECT_HANDLE_H
