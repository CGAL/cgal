// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_GENERIC_HANDLE_MAP_H
#define CGAL_GENERIC_HANDLE_MAP_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/Unique_hash_map.h>

namespace CGAL {

struct Identity_hash_function {
    std::size_t operator() (std::size_t h) const {
        return h;
    }
};

template <class I>
class Generic_handle_map : public
  Unique_hash_map<std::size_t,I,Identity_hash_function>
{ typedef Unique_hash_map<std::size_t,I,Identity_hash_function> Base;
public:
  Generic_handle_map() : Base() {}
  Generic_handle_map(I i) : Base(i) {}

  template <class H>
  std::size_t hash_function(H h) const
  { return std::size_t((void*)&*h)/sizeof(h); }

  template <class H>
  const I& operator[](H h) const
  { return Base::operator[](hash_function(h)); }

  template <class H>
  I& operator[](H h)
  { return Base::operator[](hash_function(h)); }

};

} //namespace CGAL
#endif //CGAL_GENERIC_HANDLE_MAP_H
