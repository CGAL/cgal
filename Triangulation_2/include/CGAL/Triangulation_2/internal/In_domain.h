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
#include <type_traits>

namespace CGAL {

namespace internal {

template<class T>
class Has_member_is_in_domain
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static char f(check<bool(C::*)(void) const, &C::is_in_domain>*);
  template<class C>
  static char f(check<bool(C::*)(void), &C::is_in_domain>*);

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(nullptr)) == sizeof(char));
};

template <typename F, typename FH>
inline
std::enable_if_t<Has_member_is_in_domain<F>::value, bool>
get_in_domain_impl(FH fh)
{
  return fh->is_in_domain();
}

template <typename F, typename FH>
inline
std::enable_if_t<!Has_member_is_in_domain<F>::value, bool>
get_in_domain_impl(FH )
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
      return get_in_domain_impl<Face, key_type>(k);
  }

  friend void put(In_domain, const key_type& k, const value_type& v)
  {
    k->set_in_domain(v);
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_TRIANGULATION_2_IN_DOMAIN_H
