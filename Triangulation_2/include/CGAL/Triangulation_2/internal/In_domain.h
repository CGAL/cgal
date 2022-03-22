// Copyright(c) 2022  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_TRIANGULATION_2_IN_DOMAIN_H
#define CGAL_TRIANGULATION_2_IN_DOMAIN_H

#include <CGAL/license/Triangulation_2.h>

#include <boost/property_map/property_map.hpp>

namespace CGAL {

namespace internal {

  template<class T>
class Has_member_in_domain
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static char f(check<void(C::*)(void), &C::in_domain>*);

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

  template <typename F, typename FH>
inline
typename boost::enable_if<Has_member_in_domain<F>, bool>::type
get_in_domain_impl(FH fh)
{
  return fh->in_domain();
}

template <typename F, typename FH>
inline
typename boost::disable_if<Has_member_in_domain<F>, bool>::type
get_in_domain_impl(FH fh)
{
  return false;
}



template <typename CDT>
struct In_domain {

  typedef typename CDT::Face Face;
  typedef typename CDT::Face_handle key_type;
  typedef bool value_type;
  typedef bool reference;
  typedef boost::read_write_property_map_tag category;

  friend bool get(In_domain, const key_type& k)
  {
      return get_in_domain_impl<Face, Face_handle>(k);
  }

  friend void put(In_domain, const key_type& k, const value_type& v)
  {
    k->set_in_domain(v);
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_TRIANGULATION_2_IN_DOMAIN_H
