// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_BOX_LIMITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_LIMITS_H

#include <CGAL/basic.h>
#include <CGAL/known_bit_size_integers.h>
#include <CGAL/long_long.h>
#include <cfloat>
#include <climits>

CGAL_BEGIN_NAMESPACE

namespace Box_intersection_d {


template<class T>
struct box_limits {};

template<>
struct box_limits<int> {
    static int inf() { return INT_MIN; }
    static int sup() { return INT_MAX; }
};

template<>
struct box_limits<unsigned int> {
    static unsigned int inf() { return 0; }
    static unsigned int sup() { return UINT_MAX; }
};

template<>
struct box_limits<float> {
    static float inf() { return -sup(); }
    static float sup()
    {
        const UInteger32 i = 0x7f800000;
        return reinterpret_cast<const float&>(i);
    }
};

template<>
struct box_limits<double> {
    static double inf() { return -sup(); }
    static float sup()
    {
        const UInteger64 i = 0x7FF0000000000000ull;
        return reinterpret_cast<const float&>(i);
    }
};

} // end namespace Box_intersection_d


CGAL_END_NAMESPACE


#endif
