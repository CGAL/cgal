 
// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

// NOT INTRODUCED YET 

#ifndef CGAL_POLYNOMIAL_MODULAR_TRAITS_TRAITS_H
#define CGAL_POLYNOMIAL_MODULAR_TRAITS_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Modular_traits.h>

CGAL_BEGIN_NAMESPACE

/*! \ingroup NiX_Polynomial
 *  \ingroup NiX_Modular_traits_spec
 *  \brief Specialization of Modular_traits for NiX::Polynomial.
 * 
 *  NiX::Modular_traits::Modular_image maps the coefficients of a polynomial
 *  to their Modular_image and returns the resulting polynomial.  
 */
template< class COEFF >
class Modular_traits< Polynomial<COEFF> > {
    
private:
    typedef Modular_traits<COEFF> Mtr;
public:
    typedef Polynomial<COEFF> NT;
    typedef Modular_traits<NT> Self;
    typedef typename Mtr::Is_modularizable Is_modularizable;
    typedef Polynomial<typename Mtr::Modular_NT> Modular_NT;
    
    struct Modular_image{
        Modular_NT operator()(const NT& p){ 
            std::vector<typename Mtr::Modular_NT> V;
            typename Mtr::Modular_image modular_image;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image(p[i]));
            return Modular_NT(V.begin(),V.end());           
        }
    };

    struct Modular_image_inv{ 
        NT operator()(const Modular_NT& p) const {  
            std::vector<COEFF> V;
            typename Mtr::Modular_image_inv modular_image_inv;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image_inv(p[i]));
            return NT(V.begin(),V.end());           
        }
    };
};


CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_MODULAR_TRAITS_TRAITS_H
