// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer 
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$

#ifndef CGAL_IPOWER_H
#define CGAL_IPOWER_H

#include <CGAL/assertions.h>

namespace CGAL {

template <typename NT>
inline
NT ipower(const NT& base, int expn) {
    // compute base^expn using square-and-multiply
    CGAL_precondition(expn >= 0);
    
    // handle trivial cases efficiently
    if (expn == 0) return NT(1);
    if (expn == 1) return base;
    
    // find the most significant non-zero bit of expn
    int e = expn, msb = 0;
    while (e >>= 1) msb++;
    
    // computing base^expn by square-and-multiply
    NT res = base;
    int b = 1<<msb;
    while (b >>= 1) { // is there another bit right of what we saw so far?
        res *= res;
        if (expn & b) res *= base;
    }
    return res;
}

template <typename NT>
inline
NT ipower(const NT& base, long expn) {
    // compute base^expn using square-and-multiply
    CGAL_precondition(expn >= 0);
    
    // handle trivial cases efficiently
    if (expn == 0) return NT(1);
    if (expn == 1) return base;
    
    // find the most significant non-zero bit of expn
    int e = expn, msb = 0;
    while (e >>= 1) msb++;
    
    // computing base^expn by square-and-multiply
    NT res = base;
    int b = 1<<msb;
    while (b >>= 1) { // is there another bit right of what we saw so far?
        res *= res;
        if (expn & b) res *= base;
    }
    return res;
}

} //namespace CGAL

#endif // CGAL_IPOWER_H
