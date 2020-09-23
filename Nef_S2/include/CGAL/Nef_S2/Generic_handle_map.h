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

struct Void_handle_hash_function {
    std::size_t operator() (void* h) const {
        return std::size_t(h);
    }
};


template <class I>
class Generic_handle_map : public
  Unique_hash_map<void*,I,Void_handle_hash_function>
{ typedef Unique_hash_map<void*,I,Void_handle_hash_function> Base;
public:
  Generic_handle_map() : Base() {}
  Generic_handle_map(I i) : Base(i) {}

  template <class H>
  const I& operator[](H h) const
  { return Base::operator[](&*h); }

  template <class H>
  I& operator[](H h)
  { return Base::operator[](&*h); }

};

} //namespace CGAL
#endif //CGAL_GENERIC_HANDLE_MAP_H
