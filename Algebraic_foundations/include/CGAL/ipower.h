// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// $URL: svn+ssh://afabri@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
// 
//
// Author(s)     : Michael Hemmer 
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$

#ifndef CGAL_IPOWER_H
#define CGAL_IPOWER_H

CGAL_BEGIN_NAMESPACE

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


CGAL_END_NAMESPACE

#endif // CGAL_IPOWER_H
