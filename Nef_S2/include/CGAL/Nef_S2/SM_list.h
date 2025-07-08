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
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SM_LIST_H
#define CGAL_SM_LIST_H

#include <CGAL/license/Nef_S2.h>

#include <CGAL/In_place_list.h>

namespace CGAL {

template < class T>
class SM_in_place_list
    : public T,
      public In_place_list_base<SM_in_place_list<T> > {
public:
    typedef SM_in_place_list<T> Self;
    SM_in_place_list() {}
    SM_in_place_list(const SM_in_place_list& other)=default;
    SM_in_place_list(const T& v)   // down cast
        : T(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((T*)this) = ((const T&)v);
        return *this;
    }
};

} //namespace CGAL

#endif // CGAL_SM_LIST_H
