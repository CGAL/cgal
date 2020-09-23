// Copyright (c) 1997-2000
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_SIMPLE_OBJECTS_H
#define CGAL_SIMPLE_OBJECTS_H

namespace CGAL {

template <class R>
struct Lt_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) == SMALLER; }
};

template <class R>
struct Le_from_compare {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { typename R::Compare_lexicographically_d cmp;
    return cmp(p1,p2) != LARGER; }
};

template <class R>
struct Eq_from_method {
  typedef typename R::Point_d Point_d;
  bool operator()(const Point_d& p1, const Point_d& p2) const
  { return p1 == p2; }
};

} //namespace CGAL

#endif // CGAL_SIMPLE_OBJECTS_H
