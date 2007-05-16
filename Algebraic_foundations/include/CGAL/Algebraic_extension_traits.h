// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================


/*! \file CGAL/Algebraic_extension_traits.h
 *  \brief Defines traits class CGAL::Algebraic_extension_traits. 
*/

#ifndef CGAL_ALGEBRAIC_NUMBER_TRAITS_H
#define CGAL_ALGEBRAIC_NUMBER_TRAITS_H 1

#include <numeric> // for std::accumulate

CGAL_BEGIN_NAMESPACE

template< class T >
class Algebraic_extension_traits {
public:
    //! \name Typedefs 
    //! the number type for which this instance has been instantiated
    typedef T Type;
    //! standard number types are not extended
    typedef CGAL::Tag_false Is_extended;
  
    //! computes the factor which normalizes a number to be integral after 
    //  multiplication
    class Normalization_factor 
        : public Unary_function<Type,Type> {
    private:
        static Type 
        normalization_factor(const Type&,Integral_domain_without_division_tag){
            return Type(1);
        }
        static Type 
        normalization_factor(const Type& a, Field_tag){
            return Type(1)/a;
        }
    public:
        //! determine normalization factor
        Type operator () (const Type& a) {
            CGAL_precondition(a != Type(0));
            typedef typename Algebraic_structure_traits<Type>::Algebraic_category
                Tag;
            return normalization_factor(a, Tag());
        }
    };
    
    class Denominator_for_algebraic_integers 
        : public Unary_function<Type,Type> {
    public: 
        //! determine normalization factor
        Type operator () (const Type&) {
            return Type(1);
        }
        
        template <class InputIterator>
        Type operator () (InputIterator, InputIterator) {
            return Type(1);
        }
    };
};

// This is the specialization for Sqrt_extension
// This has been moved here because Algebraic_extension_traits won't be part of 
// CGAL release 3.3
// TODO: move back to Sqrt_extension.h for release 3.4 ?

// fwd of Sqrt_extension
template <typename A, typename B> class Sqrt_extension;

template <class COEFF, class ROOT >
class Algebraic_extension_traits<Sqrt_extension<COEFF,ROOT> > {
// needed to 'add up' sqrt_extensions in iterator range such that all roots are
//   collected in order to keep operation time minimal all scalar coeffs are set
//   to 1 by standardise.
//   TODO .. find a better name, if you want to.
//
template <class NT_>
class Standardise {
public:
    typedef NT_ argument_type;
    typedef NT_ result_type;
    NT_ operator () (const NT_&) const {
        return NT_(1);
    }
};

template <class COEFF_, class ROOT_>
class Standardise<Sqrt_extension<COEFF_,ROOT_> > {
    Standardise<COEFF_> standardise;
public:
    typedef Sqrt_extension<COEFF_,ROOT_> NT_;
    typedef NT_ argument_type;
    typedef NT_ result_type;
    NT_ operator () (const NT_& a) const {
        if (a.a1() != COEFF_(0)){
            return NT_(standardise(a.a0()),standardise(a.a1()),a.root());
        }else{
            return NT_(standardise(a.a0()));
        }
    }
};

public:
    //! \name Typedefs
    //! the number type for which this instance has been instantiated
    typedef Sqrt_extension<COEFF,ROOT> NT;
    //! Sqrt_extension as a number type is extended
    typedef ::CGAL::Tag_true Is_extended;

    //! computes the factor which normalizes a number to be integral
    //!  after multiplication
    //!
    class Normalization_factor{
    public:
        //! argument type
        typedef NT argument_type;
        //! result type
        typedef NT result_type;
        //! determine normalization factor
        NT operator () (const NT& a) const {
            
            typename Algebraic_structure_traits<COEFF>::Is_zero is_zero;
            
            typedef Algebraic_extension_traits<COEFF> SET;
            typename  SET::Normalization_factor normalization_factor;
            CGAL_precondition(a != NT(0));

            NT result;
            if(a.is_extended() && ! is_zero(a.a1())){
                NT tmp1(a.a0(),-a.a1(),a.root());

                NT tmp2= a*tmp1;
                CGAL_postcondition(tmp2.a1()==COEFF(0));
                result = tmp1*NT(normalization_factor(tmp2.a0()));
                CGAL_postcondition(! is_zero (result.a1()));
            }else{
                result = NT(normalization_factor(a.a0()));
            }
            return result;
        }
    };

    //! returns the extension factor needed for the gcd_utcf computation
    //! for more details see ... TODO!!
    //!
    class Denominator_for_algebraic_integers {
    public:
        //! argument type
        typedef NT argument_type;
        //! result type
        typedef NT result_type;
        //! determine denominator for algebraic integers
    public:
        NT operator () (const NT& a) const {
            typedef Algebraic_extension_traits<COEFF> ANT;
            typename ANT::Denominator_for_algebraic_integers dfai;

            Standardise<COEFF> standardise;
            if (a.a1() != COEFF(0)) {
                COEFF tmp =
                    standardise(a.a0())
                    + standardise(a.a1())
                    + standardise(COEFF(a.root()));
                return  NT(COEFF(4) * COEFF(a.root()))* NT(dfai(tmp));
            } else {
                return NT(dfai(a.a0()));
            };
        }

        //! overloaded operator for computing the denominator out of an iterator
        //!  range accumulates all root expressions which appear in the range to
        //!  the most complex term and uses this term to determine the denominator
        //!  for algebraic integers
        //!
        template <class InputIterator>
        NT operator () (InputIterator begin, InputIterator end) const {
            NT a = std::accumulate(::boost::make_transform_iterator(begin,Standardise<NT>()),
                                   ::boost::make_transform_iterator(end  ,Standardise<NT>()), NT(0));
            return (*this)(a);
        }
    };
};

CGAL_END_NAMESPACE

#endif // NiX_ALGEBRAIC_NUMBER_TRAITS_H
// EOF
