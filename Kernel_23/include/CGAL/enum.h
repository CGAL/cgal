// Copyright (c) 1999
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
// Author(s)     : Stefan Schirra

#ifndef CGAL_ENUM_H
#define CGAL_ENUM_H

#include <CGAL/config.h>
#include <CGAL/Kernel/Same_uncertainty.h>
#include <CGAL/Origin.h>

// If you add/change one type here, please update Is_a_predicate.h as well.

namespace CGAL {

enum  Sign
{
    NEGATIVE = -1, ZERO = 0, POSITIVE = 1,

    // Orientation constants:
    RIGHT_TURN = -1, LEFT_TURN = 1,

    CLOCKWISE = -1, COUNTERCLOCKWISE = 1,

    COLLINEAR = 0, COPLANAR = 0, DEGENERATE = 0,

    // Oriented_side constants:
    ON_NEGATIVE_SIDE = -1, ON_ORIENTED_BOUNDARY = 0, ON_POSITIVE_SIDE = 1,

    // Comparison_result constants:
    SMALLER = -1, EQUAL = 0, LARGER = 1
};

typedef Sign Orientation;
typedef Sign Oriented_side;
typedef Sign Comparison_result;

enum  Bounded_side
      {
        ON_UNBOUNDED_SIDE = -1,
        ON_BOUNDARY,
        ON_BOUNDED_SIDE
      };

enum  Angle
      {
          OBTUSE = -1,
          RIGHT,
          ACUTE
      };


template <class T>
inline
T
opposite(const T& t)
{ return -t; }

inline
Sign
operator-(Sign o)
{ return static_cast<Sign>( - static_cast<int>(o)); }

inline
Bounded_side
opposite(Bounded_side bs)
{ return static_cast<Bounded_side>( - static_cast<int>(bs)); }

inline
Angle
opposite(Angle a)
{ return static_cast<Angle>( - static_cast<int>(a)); }

inline Sign operator* (Sign s1, Sign s2)
{
    return static_cast<Sign> (static_cast<int> (s1) * static_cast<int> (s2));
}


enum Box_parameter_space_2
     {
        LEFT_BOUNDARY = 0,
        RIGHT_BOUNDARY,
        BOTTOM_BOUNDARY,
        TOP_BOUNDARY,
        INTERIOR,
        EXTERIOR
     };

template < typename T, typename U >
inline
T enum_cast(const U& u)
{ return static_cast<T>(u); }

} //namespace CGAL

#endif // CGAL_ENUM_H
