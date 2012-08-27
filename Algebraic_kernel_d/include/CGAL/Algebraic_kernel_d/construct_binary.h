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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_BINARY_H
#define CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_BINARY_H

#include <CGAL/basic.h>
#include <CGAL/ipower.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#endif

#include <limits>

namespace CGAL {

namespace internal {

// Generic construct_binary function, using ipower 
template< class Integer >
inline void construct_binary( const Integer& e, Integer& x ) {
    CGAL_precondition( e >= 0 );
    Integer exponent(e);
    x = Integer(1);
    
    const Integer max_ipower = (exponent > Integer((std::numeric_limits<int>::max)())) ?
      CGAL::ipower( Integer(2), (std::numeric_limits<int>::max)() ) : 
                Integer(0);
    
    while( exponent > Integer((std::numeric_limits<int>::max)()) ) {
        x *= max_ipower;
        exponent -= Integer((std::numeric_limits<int>::max)());    
    }
    
    x *= CGAL::ipower( Integer(2), (int)CGAL::to_double(exponent) );    
}

template< class Integer, class Rational >
inline void construct_binary( const Integer& m, const Integer& e, Rational& x ) {
       Integer den(1), num;
       if(e>0) {
          construct_binary( e, num );
          num *= m;          
       }
       else {
          num = m;
          construct_binary( -e, den );
       }
       x = Rational(num, den);
}

// Specialization for LEDA
    
#ifdef CGAL_USE_LEDA
    
    // Constructs 2^e from an integer e. Needed in Descartes
    inline void construct_binary(const ::leda::integer& e, ::leda::integer& x) {
       typedef ::leda::integer  Integer;
       x = Integer(1) << e.to_long();
    }
    
    // Constructs m*2^e from two integers m,e. Needed in Descartes
    inline void construct_binary(const ::leda::integer& m, const ::leda::integer& e,
                  ::leda::rational& x) {
                    
       typedef ::leda::integer  Integer;
       typedef ::leda::rational Rational;
    
       Integer den(1);
       Integer num(m);
       if(e>0) {
          num <<= e.to_long();
       }
       else {
          den <<= (-e).to_long();
       }
       x = Rational(num, den);
    }

#endif // CGAL_USE_LEDA

// Specialization for CORE

#ifdef CGAL_USE_CORE
    
    // Constructs 2^e from an integer e. Needed in Descartes
    inline void construct_binary(const ::CORE::BigInt& e, ::CORE::BigInt& x) {
       typedef ::CORE::BigInt Integer;
       x = Integer(1) << ::CORE::ulongValue(e);
    }
    
    // Constructs m*2^e from two integers m,e. Needed in Descardes
    inline void construct_binary(const ::CORE::BigInt& m, const ::CORE::BigInt& e, 
                  ::CORE::BigRat& x) {
       typedef ::CORE::BigInt Integer;
       typedef ::CORE::BigRat Rational;
      
       Integer den(1);
       Integer num(m);
       if(e>0) {
          num <<= ::CORE::ulongValue(e);
       }
       else {
          den <<= ::CORE::ulongValue(-e);
       }
       x = Rational(num, den);
    }
    
#endif // CGAL_USE_CORE

} // namespace internal

} //namespace CGAL


#endif // CGAL_ALGEBRAIC_KERNEL_D_CONSTRUCT_BINARY_H
