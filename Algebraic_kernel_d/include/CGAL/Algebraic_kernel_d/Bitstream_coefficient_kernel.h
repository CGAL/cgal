// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ==========================================================================
#ifndef CGAL_BITSTREAM_COEFFICIENT_KERNEL_H
#define CGAL_BITSTREAM_COEFFICIENT_KERNEL_H 1

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/convert_to_bfi.h>

namespace CGAL {

namespace internal {

template <typename Coefficient_> struct Bitstream_coefficient_kernel {

    typedef Coefficient_ Coefficient;
    
    typedef typename 
        CGAL::Get_arithmetic_kernel<Coefficient_>::Arithmetic_kernel  
    Arithmetic_kernel;
    
    typedef typename Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;
    typedef typename Arithmetic_kernel::Integer Integer;
    typedef typename Arithmetic_kernel::Rational Bound;


    typedef typename CGAL::Algebraic_structure_traits<Coefficient>
        ::Is_zero  Is_zero;
    
    Is_zero is_zero_object() const {
        return Is_zero();
    }

    struct Convert_to_bfi : public std::unary_function
        <Coefficient,Bigfloat_interval> {
        
        Bigfloat_interval operator() (Coefficient c) const {
            return CGAL::convert_to_bfi(c);
        }
    };
    
    Convert_to_bfi convert_to_bfi_object() const {
        return Convert_to_bfi();
    }
        
    
}; // of class Bitstream_coefficient_kernel

} // namespace internal

} //namespace CGAL

#endif // CGAL_BITSTREAM_COEFFICIENT_KERNEL_H
