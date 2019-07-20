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

namespace CGAL {

/*! \ingroup CGAL_Polynomial
 *  \ingroup CGAL_Modular_traits_spec
 *  \brief Specialization of Modular_traits for CGAL::Polynomial.
 * 
 *  CGAL::Modular_traits::Modular_image maps the coefficients of a polynomial
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
    typedef Polynomial<typename Mtr::Residue_type> Residue_type;
    
    struct Modular_image{
        Residue_type operator()(const NT& p){ 
            std::vector<typename Mtr::Residue_type> V;
            typename Mtr::Modular_image modular_image;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image(p[i]));
            return Residue_type(V.begin(),V.end());           
        }
    };

    struct Modular_image_representative{ 
        NT operator()(const Residue_type& p) const {  
            std::vector<COEFF> V;
            typename Mtr::Modular_image_representative modular_image_representative;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image_representative(p[i]));
            return NT(V.begin(),V.end());           
        }
    };
};


} //namespace CGAL
#endif // CGAL_POLYNOMIAL_MODULAR_TRAITS_TRAITS_H
