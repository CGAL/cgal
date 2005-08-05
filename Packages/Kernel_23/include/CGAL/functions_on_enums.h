// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri
//                 Stefan Schirra

#ifndef CGAL_FUNCTIONS_ON_ENUMS_H
#define CGAL_FUNCTIONS_ON_ENUMS_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template <class T>
inline
T
opposite(const T& t)
{ return -t; }

inline
Sign
opposite(Sign o)
{ return static_cast<Sign>( - static_cast<int>(o)); }

inline
Comparison_result
opposite(Comparison_result o)
{ return static_cast<Comparison_result>( - static_cast<int>(o)); }

inline
Oriented_side
opposite(Oriented_side os)
{ return static_cast<Oriented_side>( - static_cast<int>(os)); }

inline
Bounded_side
opposite(Bounded_side bs)
{ return static_cast<Bounded_side>( - static_cast<int>(bs)); }

inline
Angle
opposite(Angle a)
{ return static_cast<Angle>( - static_cast<int>(a)); }

template < typename T, typename U >
inline
T enum_cast(const U& u)
{ return static_cast<T>(u); }

CGAL_END_NAMESPACE

#endif // CGAL_FUNCTIONS_ON_ENUMS_H
