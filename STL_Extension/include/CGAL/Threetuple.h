// Copyright (c) 1999,2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_THREETUPLE_H
#define CGAL_THREETUPLE_H

#include <CGAL/config.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Threetuple.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS "Please use `std::array` instead."
#include <CGAL/internal/deprecation_warning.h>

#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {

template < class T >
struct CGAL_DEPRECATED Threetuple
{
  typedef T value_type;

  T  e0, e1, e2;

  Threetuple()
  {}

  Threetuple(const T & a0, const T & a1, const T & a2)
    : e0(a0), e1(a1), e2(a2)
  {}
};

} //namespace CGAL

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_THREETUPLE_H
