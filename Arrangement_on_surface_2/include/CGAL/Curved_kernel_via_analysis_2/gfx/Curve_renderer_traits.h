// Copyright (c) 2004-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_CKVA_CURVE_RENDERER_TRAITS_H
#define CGAL_CKVA_CURVE_RENDERER_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/function_objects.h>

/*! \file CGAL/Curved_kernel_via_analysis_2/gfx/Curve_renderer_traits.h
 * \brief
 * defines class Curve_renderer_traits.
 * 
 * provides specialisations of Curve_renderer_traits for different number
 * types compatible with the curve renderer
*/
 
namespace CGAL {

// transformation routine
namespace internal {

//! this exception is thrown whenever the precision of used number
//! type is not sufficient
class Insufficient_rasterize_precision_exception
{  };

#ifndef CGAL_MAX_COEFF_TRANSFORM
#define CGAL_MAX_COEFF_TRANSFORM

#ifndef CGAL_ABS
#define CGAL_ABS(x) ((x) < 0 ? -(x): (x))
#endif

#ifndef CGAL_SGN
#define CGAL_SGN(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
#endif

template <class NT>
struct Max_coeff
{
    template <class X>
    NT operator()(const CGAL::Polynomial<X>& p) const
    {
        typename CGAL::Polynomial<X>::const_iterator it = p.begin();
        Max_coeff<NT> max_coeff;
        NT max(max_coeff(*it));
        while(++it != p.end()) {
            NT tmp(max_coeff(*it));
            if(max < CGAL_ABS(tmp))
                max = CGAL_ABS(tmp);
        }
        return max;
    }
    NT operator()(const NT& x) const
    { return CGAL_ABS(x); }
};

/*!\brief
 * divides an input value by a contant
 *
 * provided that there is a coercion between \c Input and \c Result types
 */
template <class Result, class Input>
struct Reduce_by {

    typedef Input argument_type;
    typedef Result result_type;

    Reduce_by(const Input& denom_) :
        denom(cast(denom_)) {
    }
    
    Result operator()(const Input& x) {
        return (cast(x)/denom);
    }

    typename CGAL::Coercion_traits<Input, Result>::Cast cast;
    Result denom;
};

/*!\brief
 * transforms bivaritate polynomial of type \c InputPoly_2 to
 * \c OutputPoly_2 by recursively applying operation \c Op to all of its
 * coefficients
 *
 * <tt>Op: InputPoly_2::Inntermost_coefficient_type ->
 *             OutputPoly_2::Inntermost_coefficient_type</tt>
 */
template <class OutputPoly_2, class InputPoly_2, class Op>
struct Transform {

    typedef InputPoly_2  first_argument_type;
    typedef Op           second_argument_type;
    typedef OutputPoly_2 result_type;

    template <class X>
    OutputPoly_2 operator()(const CGAL::Polynomial<X>& p, Op op = Op()) const {

        Transform<typename OutputPoly_2::NT, typename InputPoly_2::NT, Op> tr;
        return OutputPoly_2(
            ::boost::make_transform_iterator(p.begin(), std::bind2nd(tr, op)),
            ::boost::make_transform_iterator(p.end(), std::bind2nd(tr, op)));
    }

    OutputPoly_2 operator()(
        const typename Innermost_coefficient_type<InputPoly_2>::Type& x, Op op)
        const {
        return static_cast<OutputPoly_2>(op(x));
    }
};

#endif // CGAL_MAX_COEFF_TRANSFORM

/*!\brief
 * class template \c Curve_renderer_traits
 *
 * this traits class prodives various number type conversions for the
 * curve renderer
 */              

template <class Coeff_, class Integer_, class Rational_>
struct Curve_renderer_traits_base 
{ 
    //! type of innermost polynomial coefficients
    typedef Coeff_ Coeff;
    
    //! an integer number type
    typedef Integer_ Integer; 

    typedef Rational_ Rational;

    //! conversion from rational to floating-point
    template <class Float>
    struct Rat_to_float {
        typedef Float result_type;

        template <class X, class Y,class ACDE_TAG,class FP_TAG>
        Float operator()(const Sqrt_extension<X, Y, ACDE_TAG, FP_TAG>& x) const { 
            typename CGAL::Coercion_traits<Sqrt_extension<X, Y, ACDE_TAG, FP_TAG>, Float>::Cast
                cast;        
            return cast(x); 
        }

        Float operator()(const Rational& x) const { 
            return static_cast<Float>(x); 
        }
    };
   
    //! provided for convenience when there exists an implicit coercion
    //! between number types
    template <class To>
    struct Implicit_coercion {
        typedef To result_type;
        
        template <class From> 
        To operator()(const From& x) const {
            return static_cast<To>(x);
        }
    };
   
    struct Float_to_int {
        typedef int result_type;
        
        template <class Float>
        int operator()(const Float& x) const
        { return static_cast<int>(std::floor(CGAL::to_double(x))); }
            //return static_cast<int>(std::floor(x)); }
    };

    /*!\brief
     * conversion from bivariate polynomial over extended number type to
     * polynomial with coefficients of type \c Coeff
     *
     * valid instantiations of \c Extended number type are \c Rational ,
     * \c FieldWithSqrt , etc. Provided that there is a type coercion between
     * \c Extended and \c Coeff
     */
    struct Convert_poly { 
        typedef CGAL::Polynomial<CGAL::Polynomial<Coeff> > result_type;
        
        template <class Extended>
        inline result_type operator()(const 
            CGAL::Polynomial<CGAL::Polynomial<Extended> >& poly) const {
            
    //!@todo: use Rat_to_float functor instead of coercion traits ?
    //! need some sort of polymorphism..

            //std::cerr << "calling transform..\n";
            typedef typename CGAL::Coercion_traits<Extended, Coeff>::Cast
                Cast;
            Transform<result_type,
                    CGAL::Polynomial<CGAL::Polynomial<Extended> >, Cast>
                transform;
            return transform(poly);
        }
    };
    
    //! converts polynomial coefficients to floating-point representation
    //! \c error_bounds is set if the result is not reliable
    struct Extract_eval {
        typedef Coeff argument_type;
        typedef Coeff result_type;
        
        Coeff operator()(const Coeff& x, 
            bool *error_bounds = NULL) const { 
            if(error_bounds != NULL)
                 *error_bounds = false;
            return x;
        }
    };
    
    //! computes a 32-bit hash value from floating-point argument
    struct Hash_function {
        typedef std::size_t result_type;
        struct long_long {
            long low, high;
        };
        
        template <class Float>
        std::size_t operator()(const Float& key) const {
            const long_long *hk = reinterpret_cast<const long_long *>(&key);
            return (hk->low ^ hk->high);
        }
    };
    
    //! makes result exact after inexact operation such as div, sqrt, etc.
    //! (required for multi-precision arithmetic)
    struct Make_exact {
        typedef void result_type;
        
        template <class Float>
        void operator()(const Float& x) const
        { }
    };
    
    //! compares a given quantity with the precision limit a floating-point
    //! number type, returns \c true if this limit is exceeded
    struct Precision_limit {
        typedef bool result_type;
        
        template <class Float>
        bool operator()(const Float& x) const
        { return false;/*(CGAL_ABS(x) <= 1e-16)*/; }
    };
    
    //! maximum subdivision level for the curve renderer by exceeding which
    //! an exception of type \c Insufficient_rasterize_precision_exception
    //! will be thrown, this is also limited by \c Integer number type, since
    //! the integer must be able to store numbers up to 2^MAX_SUBDIVISION_LEVEL
    static const unsigned MAX_SUBDIVISION_LEVEL = 12;   
};

template <class Float, class Rational>
struct Curve_renderer_traits 
{  };

#ifdef CGAL_USE_CORE
//! Specialization for \c CGAL::Interval_nt<true>
template <>
struct Curve_renderer_traits<CGAL::Interval_nt<true>, CORE::BigRat > :
        public Curve_renderer_traits_base<CGAL::Interval_nt<true>, int,
            CORE::BigRat> {
 
    typedef double Float;

    struct Rat_to_float {
        typedef Float result_type;
                
        template <class Extended>
        Float operator()(const Extended& x) const {
            return CGAL::to_double(x); 
        }
    };

    typedef Implicit_coercion<double> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;

        Integer operator()(const Rational& x) const {
            return static_cast<int>(std::floor(CGAL::to_double(x)));
        }
    };

    struct Extract_eval {
        typedef Coeff argument_type;
        typedef Float result_type;
        
        Float operator()(const Coeff& x, 
                    bool *error_bounds = NULL) const { 
            bool err_bnd;
//             err_bnd = (CGAL_ABS(l) < 1E-15 || CGAL_ABS(u) < 1E-15) ||
//                 ((l <= 0 && u >= 0));
            Float l = x.inf(), u = x.sup(), mid = (l+u)/2;
            err_bnd = ((l < 0 && u > 0)||(l == 0 && u == 0));
            if(error_bounds != NULL)
                *error_bounds = err_bnd;
//! ATTENTION!!! if smth is screwed up try to uncomment the line below
//! this is very crucial for performance & stability
            if(err_bnd)  // &&  ABS(mid) < 1E-15)
                return 0; 
//! ATTENTION!!!
            return mid;
        }
    };

    //! compares a given quantity with the precision limit a floating-point
    //! number type, returns \c true if this limit is exceeded
    struct Precision_limit {
        typedef bool result_type;
        
        template <class Float>
        bool operator()(const Float& x) const
        { return (CGAL_ABS(x) <= 1e-16); }
    };

    static const unsigned MAX_SUBDIVISION_LEVEL = 12;
};

//! Specialization for \c CORE::BigFloat
template <>
struct Curve_renderer_traits<CORE::BigFloat, class CORE::BigRat> 
         : public Curve_renderer_traits_base<CORE::BigFloat, CORE::BigInt,
                CORE::BigRat> {

    typedef CORE::BigFloat Float;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.BigIntValue(); 
        }
    };

    typedef Rat_to_float<Float> Rat_to_float;
    
    struct Float_to_rat {
        typedef Float argument_type;
        typedef Rational result_type;
        
        Rational operator()(const Float& x) const
        { return x.BigRatValue(); }
    };
    
    struct Hash_function {
        typedef Float argument_type;
        typedef std::size_t result_type;
        
        inline result_type operator()(const Float& key) const {
            const CORE::BigFloatRep& rep = key.getRep();
            std::size_t ret = reinterpret_cast<std::size_t>(&rep);
            return ret;
        }
    };
    
    struct Make_exact {
        typedef Float argument_type;
        typedef void result_type;
        
        inline void operator()(Float& x) const
        { x.makeExact(); }
    };

    static const unsigned MAX_SUBDIVISION_LEVEL = 12;
};

//! Specialization for \c CORE::BigRat
template <>
struct Curve_renderer_traits<CORE::BigRat, CORE::BigRat> : 
    public Curve_renderer_traits_base<CORE::BigRat, CORE::BigInt,
        CORE::BigRat> {

    typedef CORE::BigRat Float;

    typedef Rat_to_float<Float> Rat_to_float;

    typedef Implicit_coercion<Float> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return x.BigIntValue(); 
        }
    };
    
    struct Hash_function {
        typedef Float argument_type;
        typedef std::size_t result_type;
        
        inline result_type operator()(const Float& key) const {
            const CORE::BigRatRep& rep = key.getRep();
            std::size_t ret = reinterpret_cast<std::size_t>(&rep);
            return ret;
        }
    };
};
#endif // CGAL_USE_CORE

#ifdef CGAL_USE_LEDA
template <>
struct Curve_renderer_traits<CGAL::Interval_nt<true>, leda::rational > :
        public Curve_renderer_traits_base<CGAL::Interval_nt<true>, int,
            leda::rational> {
 
    typedef double Float;

    struct Rat_to_float {
        typedef Float result_type;
                
        template <class Extended>
        Float operator()(const Extended& x) const {
            return CGAL::to_double(x); 
        }
    };

    typedef Implicit_coercion<double> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;

        Integer operator()(const Rational& x) const {
            return static_cast<int>(std::floor(CGAL::to_double(x)));
            //return static_cast<int>(x.to_double());
        }
    };

    struct Extract_eval {
        typedef Coeff argument_type;
        typedef Float result_type;
        
        Float operator()(const Coeff& x, 
                    bool *error_bounds = NULL) const { 
            bool err_bnd;
//             err_bnd = (CGAL_ABS(l) < 1E-15 || CGAL_ABS(u) < 1E-15) ||
//                 ((l <= 0 && u >= 0));
            Float l = x.inf(), u = x.sup(), mid = (l+u)/2;
            err_bnd = ((l < 0 && u > 0)||(l == 0 && u == 0));
            if(error_bounds != NULL)
                *error_bounds = err_bnd;
//! ATTENTION!!! if smth is screwed up try to uncomment the line below
//! this is very crucial for performance & stability
            if(err_bnd)  // &&  ABS(mid) < 1E-15)
                return 0; 
//! ATTENTION!!!
            return mid;
        }
    };

    //! compares a given quantity with the precision limit a floating-point
    //! number type, returns \c true if this limit is exceeded
    struct Precision_limit {
        typedef bool result_type;
        
        template <class Float>
        bool operator()(const Float& x) const
        { return (CGAL_ABS(x) <= 1e-16); }
    };
};

//! Specialization for \c leda::bigfloat
template <>
struct Curve_renderer_traits<leda::bigfloat, class leda::rational> 
         : public Curve_renderer_traits_base<leda::bigfloat, leda::integer,
                leda::rational> {

    typedef leda::bigfloat Float;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return leda::floor(x); 
        }
    };

    struct Rat_to_float {
        typedef Float result_type;

        template <class X, class Y,class ACDE_TAG,class FP_TAG>
        Float operator()(const Sqrt_extension<X, Y, ACDE_TAG, FP_TAG>& x) const { 
            typename CGAL::Coercion_traits<Sqrt_extension<X, Y, ACDE_TAG, FP_TAG>, Float>::Cast
                cast;        
            return cast(x); 
        }
        // no implicit coercion between leda rational and floats, therefore
        // decompose rational to compute the result
        Float operator()(const Rational& x) const { 
            return (static_cast<Float>(x.numerator()) / 
                    static_cast<Float>(x.denominator())); 
        }
    };
    
    struct Float_to_rat {
        typedef Float argument_type;
        typedef Rational result_type;
        
        Rational operator()(const Float& x) const
        { return x.to_rational(); }
    };
    
    struct Convert_poly { 
        typedef CGAL::Polynomial<CGAL::Polynomial<Coeff> > result_type;
        
        template <class Extended>
        inline result_type operator()(const 
            CGAL::Polynomial<CGAL::Polynomial<Extended> >& poly) const {
            
            Transform<result_type,
                CGAL::Polynomial<CGAL::Polynomial<Extended> >,
                    Rat_to_float> transform;
            return transform(poly);
        }
    };
    
    struct Hash_function {
        typedef Float argument_type;
        typedef std::size_t result_type;
        
        inline result_type operator()(const Float& key) const {
           return static_cast<std::size_t>(
                    key.get_significant_length());
        }
    };
};

//! Specialization for \c leda::rational
template <>
struct Curve_renderer_traits<leda::rational, leda::rational> : 
    public Curve_renderer_traits_base<leda::rational, leda::integer,
        leda::rational> {

    typedef leda::rational Float;

    typedef Rat_to_float<Float> Rat_to_float;

    typedef Implicit_coercion<Float> Float_to_rat;

    struct Rat_to_integer {
        typedef Rational argument_type;
        typedef Integer result_type;
        
        Integer operator()(const Rational& x) const { 
            return leda::floor(x);  
        }
    };
    
    struct Hash_function {
        typedef Float argument_type;
        typedef std::size_t result_type;
        
        inline result_type operator()(const Float& key) const {
            std::size_t ret = reinterpret_cast<std::size_t>(
                key.numerator().word_vector());
            return ret;
        }
    };

    struct Make_exact {
        typedef Float argument_type;
        typedef void result_type;
        
        inline void operator()(Float& x) const
        { x.normalize(); }
    };
};

#endif // CGAL_USE_LEDA

} // namespace internal

} //namespace CGAL

#endif // CGAL_CKVA_CURVE_RENDERER_TRAITS_H
