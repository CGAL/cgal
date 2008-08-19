// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_ACK_SHEAR_H
#define CGAL_ACK_SHEAR_H 1

#include <CGAL/basic.h>
#include <CGAL/Polynomial_traits_d.h>

#include <utility>
#include <vector>
#include <functional>
#include <iterator>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*! \ingroup NiX_bivariate_polynomial_hacks
 *  \brief Computes the polynomial f(x+sy,y)
 */ 
template<class NT>
CGAL::Polynomial<CGAL::Polynomial<NT> > 
shear(const CGAL::Polynomial<CGAL::Polynomial<NT> >& f,NT s) {
    typedef CGAL::Polynomial<NT> Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;

    Poly_1 x(NT(0),NT(1));
    Poly_1 zero(NT(0));
    Poly_1 one(NT(1));
    Poly_2 for_x(x,Poly_1(NT(s)));
    Poly_2 for_y(zero,one);

    std::vector<Poly_2> coeffs;
    coeffs.push_back(for_x);
    coeffs.push_back(for_y);

    return typename CGAL::Polynomial_traits_d<Poly_2>::Substitute()
        (f,coeffs.begin(), coeffs.end());
    
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // NiX_BIVARIATE_POLYNOMIAL_HACKS_H
// EOF
