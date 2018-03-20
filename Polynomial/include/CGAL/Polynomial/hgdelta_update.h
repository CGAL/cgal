// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Hemmer 
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_HGDELTA_UPDATE_H
#define CGAL_POLYNOMIAL_HGDELTA_UPDATE_H

namespace CGAL {

// This subroutine has been retained here for use in both new files.
namespace internal {
    template <class NT> inline
    void hgdelta_update(NT& h, const NT& g, int delta) {
        typename Algebraic_structure_traits<NT>::Integral_division idiv;
    
        // compute h = h^(1-delta) * g^delta
        switch (delta) {
        case 0:
            // h = h;
            break;
        case 1:
            h = g;
            break;
        default:
            h = idiv(CGAL::ipower(g, delta), CGAL::ipower(h, delta-1));
            break;
        }
    }
} // namespace internal

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_HGDELTA_UPDATE_H
