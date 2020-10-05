// Copyright (c) 2000
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_RANDOM_POLYGON_TRAITS_2_H
#define CGAL_RANDOM_POLYGON_TRAITS_2_H

namespace CGAL {

//-----------------------------------------------------------------------//
//                          Random_polygon_traits_2
//-----------------------------------------------------------------------//

template <class R_>
class Random_polygon_traits_2
{
  public:
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Less_xy_2         Less_xy_2;
    typedef typename R::Orientation_2     Orientation_2;

    Less_xy_2
    less_xy_2_object() const
    { return Less_xy_2(); }

    Orientation_2
    orientation_2_object() const
    { return Orientation_2(); }
};

} //namespace CGAL

#endif // CGAL_RANDOM_POLYGON_TRAITS_2_H
