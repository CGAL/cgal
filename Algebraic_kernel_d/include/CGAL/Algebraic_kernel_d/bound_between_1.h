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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_BOUND_BETWEEN_1_H
#define CGAL_BOUND_BETWEEN_1_H 1

#include <iterator>

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/enums.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Bigfloat_interval_traits.h>
#include <boost/numeric/interval.hpp>

#include <CGAL/Coercion_traits.h>

namespace CGAL {

namespace internal {

/*! \brief Tries to find a SIMPLE rational q with a<q<b.
 *
 * In this context, simple means that the denominator of <tt>q</tt>
 * is a power of two, and is not too big. There is no guarantee to find
 * the rational value between <tt>a</tt> and <tt>b</tt> of minimal
 * bit size.
 */
template<typename Algebraic_real>
typename Algebraic_real::Rational
simple_bound_between(const Algebraic_real& a,
                        const Algebraic_real&b) {

    //srb.start();
    typedef typename Algebraic_real::Rational Rational;
    typename CGAL::Fraction_traits<Rational>::Compose compose;
    typedef typename 
        CGAL::Get_arithmetic_kernel<Rational>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval Bigfloat_interval;
    typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval>
        ::Bound Bigfloat;
    typedef typename AK::Integer Integer;

    long old_prec = CGAL::get_precision(Bigfloat_interval());

    CGAL_assertion(a!=b);
    if(a>b) {
      return simple_bound_between(b,a);
    }
      
    //std::cout << "Intermediate1: " << CGAL::to_double(a) << " " << CGAL::to_double(b) << std::endl;
    /*
     * First, refine a and b until their isolating intervals are disjoint
     * Therefore, the bigger interval is refined in each substep
     */
    //srb_a.start();
    if(a.high() >= b.low()) {
        Rational size_a=a.high()-a.low(),
            size_b=b.high() - b.low();
        while(a.high() >= b.low()) {
            if(size_a < size_b) {
                b.refine();
                size_b=b.high() - b.low();
            } else {
                a.refine();
                size_a=a.high()-a.low();
            }
        }
    }
    //srb_a.stop();

    //srb_b.start();
    Bigfloat x=CGAL::upper(CGAL::convert_to_bfi(a.high()));
    Bigfloat y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
    
    if(x>=y) {
        Rational size_a=a.high() - a.low(),
            size_b=b.high() - b.low(),
            size_max = size_a>size_b ? size_a : size_b,
            size_int = b.low()-a.high();
        while(x>=y) {
            //std::cout << "x and y: " <<  x << " and " << y << std::endl;
            //std::cout << "sizes: " << CGAL::to_double(size_int) << " " << CGAL::to_double(size_max) << std::endl;
            if(size_int>size_max) {
                CGAL::set_precision(Bigfloat_interval(),
                                   2*CGAL::get_precision(Bigfloat_interval()));
                x=CGAL::upper(CGAL::convert_to_bfi(a.high()));
                y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
            } else {
                if(size_a < size_b) {
                    b.refine();
                    size_b=b.high() - b.low();
                    y=CGAL::lower(CGAL::convert_to_bfi(b.low()));
                } else {
                    a.refine();
                    size_a=a.high()-a.low();    
                    x=CGAL::upper(CGAL::convert_to_bfi(a.high()));
                }
                size_max = size_a>size_b ? size_a : size_b;
                size_int = b.low()-a.high();
            }
        }
    }
    CGAL_assertion(x<y); 

    //srb_b.stop();
    //std::cout << "Intermediate2: " << x << " " << y << std::endl;
    typename CGAL::internal::Float_traits<Bigfloat>::Get_mantissa mantissa;
    typename CGAL::internal::Float_traits<Bigfloat>::Get_exponent exponent;

    // std::cout << CGAL::to_double(x) << " < " << CGAL::to_double(y) << std::endl;

    Integer x_m = mantissa(x), y_m=mantissa(y);
    long x_e = exponent(x), y_e = exponent(y);
    //std::cout << "Floats1: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;
    
    
    if (((x_m > 0) && (y_m < 0)) || ((x_m < 0) && (y_m > 0))) {
        //srb.stop();
        return Rational(0);
    }
    bool negative=false;
    if(x_m<=0 && y_m <=0) {
        x_m=-x_m;
        y_m=-y_m;
        std::swap(x_m,y_m);
        std::swap(x_e,y_e);
        negative=true;
    }
    // Now, we have that (x_m,x_e) represents a number smaller than (y_m,y_e)
    //srb_c.start();
    //std::cout << "Floats2: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // As long as the mantissa is even, simplify
    while(x_m != 0 && (x_m & 1)==0 ) {
        x_m=x_m >> 1;
        x_e++;
    }
    while(y_m != 0 && (y_m & 1)==0 ) {
        y_m=y_m >> 1;
        y_e++;
    }
    //srb_c.stop();
    //std::cout << "Floats3: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // Bring both numbers to a common exponent
    //srb_d.start();
    long min_e = x_e < y_e ? x_e : y_e;
    while(x_e > min_e) {
        x_m=x_m << 1;
        x_e--;
    }
    while(y_e > min_e) {
        y_m=y_m << 1;
        y_e--;
    }
    //srb_d.stop();
    CGAL_assertion(y_e==x_e && x_e==min_e);
    CGAL_assertion(x_m < y_m); 
    //std::cout << "Floats4: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;

    // Avoid mantissas to have difference one
    if(y_m-x_m==Integer(1)) {
        x_m=x_m << 1;
        y_m=y_m << 1;
        x_e--;
        y_e--;
        min_e--;
    }
    //std::cout << "Floats5: " << x_m << " " << x_e << " and " << y_m << " " << y_e << std::endl;
    Integer final_mantissa(0);
    //srb_e.start();
    long x_log = x_m==Integer(0) ? -1 : CGAL::internal::floor_log2_abs(x_m),
        y_log = y_m==Integer(0) ? -1 : CGAL::internal::floor_log2_abs(y_m),
        old_log = y_log;
    //std::cout << x_log << " < " << y_log << std::endl;  
    while(x_log==y_log) {
        //std::cout << "here" << std::endl;
        while(old_log > y_log) {
            final_mantissa = final_mantissa << 1;
            old_log--;
        }
        CGAL_assertion((x_m & ((Integer(1) << x_log) - 1)) == x_m - CGAL::ipower(Integer(2),x_log));
        x_m = x_m & ((Integer(1) << x_log) - 1); // x_m - CGAL::ipower(Integer(2),x_log);
        y_m = y_m & ((Integer(1) << y_log) - 1); // y_m - CGAL::ipower(Integer(2),y_log);

        final_mantissa++;
        old_log=y_log;
        x_log = x_m==0 ? -1 : CGAL::internal::floor_log2_abs(x_m);
        y_log = y_m==0 ? -1 : CGAL::internal::floor_log2_abs(y_m);
    }
    //srb_e.stop();
    // Now, x_log != y_log, in fact, y_log is greater
    CGAL_assertion(x_log<y_log);
    //srb_f.start();
    while(old_log > y_log) {
        final_mantissa = final_mantissa << 1;
        old_log--;
    }
    if((y_m & ((Integer(1) << y_log) - 1 ))==0) { // y_m - CGAL::ipower(Integer(2),y_log)==0) {
        // Now, the constructed value would be equal to
        while(y_log!=0 && x_log==y_log-1) {
            final_mantissa = final_mantissa << 1;
            final_mantissa++;
            y_log--;
            x_m = x_m==0 ? 0 : x_m & ((Integer(1) << x_log) - 1); //x_m - CGAL::ipower(Integer(2),x_log);
            x_log = x_m==0 ? -1 : CGAL::internal::floor_log2_abs(x_m);
        }
        final_mantissa = final_mantissa << 1;
        final_mantissa++;
        y_log--;
    } else {
        final_mantissa++;
    }
    //srb_f.stop();
    min_e += y_log;
    Rational rat_between;
    //std::cout << "Min_e: " << min_e << std::endl;
    if(min_e > 0) {
        rat_between = compose(final_mantissa << min_e,
                              Integer(1));
    } else {
        rat_between = compose(final_mantissa, Integer(1) << -min_e);
    }
    if(negative) {
        rat_between = -rat_between;
    }
    //std::cout << "Result: " << a.high() << " " << rat_between << " " << b.low() << std::endl;
    CGAL_assertion(a.high() < rat_between);
    CGAL_assertion(b.low() > rat_between);
    CGAL::set_precision(Bigfloat_interval(),old_prec);
    //srb.stop();
    return rat_between;
}

	  
} // namespace internal


} //namespace CGAL

#endif // CGAL_BOUND_BETWEEN_1_H
