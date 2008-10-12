// Copyright (c) 2003-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_IOSTREAM_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_IOSTREAM_H 1

/*! \file include/CGAL/Algebraic_kernel_d/IO/Algebraic_surface_3_iostream.h
 * \brief IO operator for Algebraic_surface_3
 */

#include <iostream>
#include <fstream>

#include <CGAL/Polynomial_type_generator.h>
//#include <CGAL/Polynomial/polynomial_functions.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_surface_3.h>

#include <CGAL/Benchmark/benchmark_format.hpp>

CGAL_BEGIN_NAMESPACE

/*! \brief extracts a trivariate polynomial from a stream
 */
template < class Coefficient >
void input_polynomial_3(
        std::istream& is, 
        CGAL::Polynomial< CGAL::Polynomial< CGAL::Polynomial< Coefficient > > >& p) { 
    
    CGAL_precondition(CGAL::get_mode(is) != CGAL::IO::PRETTY);
    
    //! type of defining polynomial
    typedef typename 
        CGAL::Polynomial_type_generator< Coefficient, 3 >::Type Poly_coeff_3;
    
    Poly_coeff_3 input_poly;
    
    switch (CGAL::get_mode(is)) {
    case CGAL::IO::ASCII:
        char c;
        is >> c;
        if (c == 'P' || c == 'p') {
            is >> input_poly;
        }
        break;
    case ::CGAL::IO::BINARY:
        CGAL_error_msg("Binary input for Algebraic_surface_3 not implemented");
        break;
    default:
        CGAL_error_msg("Pretty input for Algebraic_surface_3 not allowed");
        break;
    }
    
    p = input_poly;
    
    return;
}

#if 0

template< class ArithmeticKernel >
class Benchmark_rep< CGAL::Algebraic_surface_3< ArithmeticKernel > > {
    const CGAL::Algebraic_surface_3< ArithmeticKernel >& t;
public:
    //! initialize with a const reference to \a t.
    Benchmark_rep( const CGAL::Algebraic_surface_3< ArithmeticKernel >& tt) :
        t(tt) {
    }
    
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const { 
        CGAL_error_msg("Benchmark output for Algebraic_surface_3 not implemented");
        return out << "Algebraic_surface_3(" 
#if 0
                   << bmformat( CGAL::coefficient(t.f(), 2,0,0) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 1,1,0) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 1,0,1) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,2,0) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,1,1) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,0,2) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 1,0,0) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,1,0) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,0,1) ) << ","
                   << bmformat( CGAL::coefficient(t.f(), 0,0,0) )
#endif
                   << ")";
    }
};

#endif
    
/*! \brief inserts \a quadric into output stream \a os
 *  \relates Quadric_3
 */
template < class ArithmeticKernel >
std::ostream& operator<<(
        std::ostream& os, 
        const CGAL::Algebraic_surface_3< ArithmeticKernel >& surface) {
    
    switch (CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        os << surface.f();
        break;
    case CGAL::IO::ASCII:
        // write only polynomial (P). We do not know interpolation points (I)
        os << "P ";
        os << surface.f();
        break;
    case CGAL::IO::BINARY:
        CGAL_error_msg("Binary output for Algebraic_surface_3 not implemented");
        break;
    default:
        break;
    }
    
    return os;
}


/*! \brief extracts an algebraic surface from the input 
 * stream \a is and stores it in \a surface
 */
template < class ArithmeticKernel >
std::istream& operator>>(
        std::istream& is, 
        CGAL::Algebraic_surface_3< ArithmeticKernel >& surface) { 
    
    typedef ArithmeticKernel Arithmetic_kernel;

    //! type of integer
    typedef typename Arithmetic_kernel::Integer Integer;

    //! type of defining polynomial
    typedef typename 
        CGAL::Polynomial_type_generator< Integer, 3 >::Type Poly_int_3;
    
    Poly_int_3 input_poly;
    
    CGAL::input_polynomial_3< Integer >(is, input_poly);
    
    surface = Algebraic_surface_3< ArithmeticKernel >(input_poly);
    
    return is;
}


template < class ArithmeticKernel, class OutputIterator >
bool read_file(const char *filename, OutputIterator result) {

    typedef ArithmeticKernel Arithmetic_kernel;
    
    //! type of integer
    typedef typename Arithmetic_kernel::Integer Integer;

    //! type of defining polynomial
    typedef typename
        CGAL::Polynomial_type_generator< Integer, 3 >::Type Poly_int_3;

    std::ifstream file(filename);
    
    if (!file) {
        return false;
    }
    
    std::string header;
    file >> header;

    int n = atoi(header.c_str());
    Poly_int_3 p;
    std::vector< Poly_int_3 > temp;
    for (int i = 0; i < n; i++) {
        CGAL::input_polynomial_3< Integer >(file, p);
        //*result++ = p;
        temp.push_back(p);
    }
    
    std::copy(temp.begin(), temp.end(), result);
    
    file.close();
    return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_IOSTREAM_H
// EOF
