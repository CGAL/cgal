// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Constantinos Tsirogiannis

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_CONVERTER_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_CONVERTER_H

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/NT_converter.h>

namespace CGAL {

// TODO :
// - FT converter ?

template < class Al_K1, class Al_K2,
           class RT_converter = NT_converter<typename Al_K1::RT,
                                             typename Al_K2::RT>,
           class Root_of_converter = NT_converter<typename Al_K1::Root_of_2,
                                                  typename Al_K2::Root_of_2 > >
class Algebraic_kernel_converter {
public:

        typedef typename Al_K1::RT      RT_1;
        typedef typename Al_K2::RT      RT_2;
        typedef RT_converter            RT_type_converter;
        typedef Root_of_converter       Root_of_type_converter;

        typename Al_K2::Polynomial_1_2 operator () (const typename Al_K1::Polynomial_1_2 &p) const
        {
                return typename Al_K2::Polynomial_1_2(RT_converter()(p.a()),
                                                      RT_converter()(p.b()),
                                                      RT_converter()(p.c()));
        }

        typename Al_K2::Polynomial_for_circles_2_2 operator ()
                       (const typename Al_K1::Polynomial_for_circles_2_2 &p) const
        {
                return typename Al_K2::Polynomial_for_circles_2_2(RT_converter()(p.a()),
                                                                  RT_converter()(p.b()),
                                                                  RT_converter()(p.r_sq()));
        }
};

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_CONVERTER_H
