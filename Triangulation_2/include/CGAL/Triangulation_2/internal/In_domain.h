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

template <typename CDT>
struct In_domain {

  typedef typename CDT::Face_handle key_type;
  typedef bool value_type;
  typedef bool reference;
  typedef boost::read_write_property_map_tag category;

  friend bool get(In_domain, const key_type& k)
  {
    return k->is_in_domain();
  }

  friend void put(In_domain, const key_type& k, const value_type& v)
  {
    k->set_in_domain(v);
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_TRIANGULATION_2_IN_DOMAIN_H
