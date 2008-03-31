// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_HGDELTA_UPDATE_H
#define CGAL_POLYNOMIAL_HGDELTA_UPDATE_H

CGAL_BEGIN_NAMESPACE

// This subroutine has been retained here for use in both new files.
namespace CGALi {
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
            h = idiv(ipower(g, delta), ipower(h, delta-1));
            break;
        }
    }
} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_HGDELTA_UPDATE_H
