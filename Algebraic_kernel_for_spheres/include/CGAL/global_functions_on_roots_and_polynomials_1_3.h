// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {

template < class AK >
inline
typename AK::Polynomial_1_3
construct_polynomial_1_3(const typename AK::FT& a,
                         const typename AK::FT& b,
                         const typename AK::FT& c,
                         const typename AK::FT& d)
{ return AK().construct_polynomial_1_3_object()(a, b, c, d); }

template < class AK >
inline
typename AK::Polynomials_for_line_3
construct_polynomials_for_line_3(const typename AK::FT& a1,
                                 const typename AK::FT& b1,
                                 const typename AK::FT& a2,
                                 const typename AK::FT& b2,
                                 const typename AK::FT& a3,
                                 const typename AK::FT& b3)
{ return AK().construct_polynomials_for_line_3_object()(a1, b1, a2, b2, a3, b3); }

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H
