// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_FRACTION_TRAITS_H
#define CGAL_SQRT_EXTENSION_FRACTION_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {




//################################# CGAL::Fraction_traits ##################
// Select the right alternative as Fraction_traits
// The actual Type traits is Intern::Sqrt_ext_Ftr_base_2
// The selction is done in two steps:
// 1. Inter::Sqrt_ext_Ftr_base_1 selects by the BOOL_TAG whether the COEFF type
//    Is_fraction
// 2. Intern::Sqrt_ext_Ftr_base_2 checks whether the internal type of the ROOT
//    is still implicite convertible to the new COEFF type.
//    since the ROOT type it self can not be converted.
namespace Intern{
    template <class EXT, bool> class Sqrt_ext_Ftr_base_2;
    template <class EXT, class BOOL_TAG> class Sqrt_ext_Ftr_base_1;
}

/*! \ingroup CGAL_Sqrt_extension
    \ingroup CGAL_Fraction_traits_spec
    \brief Specialisation of CGAL::Fraction_traits for CGAL::Sqrt_extension.
 *
 *  Extensions provide suitable specializations of \c CGAL::Fraction_traits.
 *  They are decomposable iff their coefficient type is.
 *  The denominator \e d of a Extension \e ext is a low common multiple
 *  (see \c CGAL::Fraction_traits::Common_factor for details) of the
 *  denominators of its coefficients.  The numerator is the Extenion
 *  \e d*ext with a fraction-free coefficient type.
 *
 *  This works for nested Sqrt_extensions, too.
 */

template <class COEFF, class ROOT, class ACDE_TAG, class FP_TAG >
class Fraction_traits< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG > >
    : public Intern::Sqrt_ext_Ftr_base_1<
    Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG >,
    typename CGAL::Fraction_traits<COEFF>::Is_fraction >
{
    // nothing new
};

namespace Intern {

// Use this if the coefficients cannot be decomposed
// into numerator and denominator
template <class NT_ >
class Sqrt_ext_Ftr_base_2< NT_, false > {
public:
    typedef NT_ NT;
    typedef ::CGAL::Tag_false Is_fraction;
    typedef ::CGAL::Null_tag Numerator_type;
    typedef ::CGAL::Null_tag Denominator_type;
    typedef ::CGAL::Null_tag Common_factor;
    typedef ::CGAL::Null_tag Decompose;
    typedef ::CGAL::Null_tag Compose;
};

template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
class Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>, true > {
private:
    typedef Fraction_traits<COEFF> CFT;
public:
  typedef Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> NT;
    typedef CGAL::Tag_true Is_fraction;
  typedef Sqrt_extension<typename CFT::Numerator_type,ROOT,ACDE_TAG,FP_TAG> Numerator_type;
    typedef typename CFT::Denominator_type Denominator_type;
    typedef typename Algebraic_structure_traits<Denominator_type>::Gcd Common_factor;

    class Decompose {
    public:
        typedef NT first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Denominator_type& third_argument_type;
        void operator () (const NT& ext,
                          Numerator_type&   num,
                          Denominator_type& den){
            typename CFT::Decompose decompose;
            typename CFT::Common_factor common_factor;
            typedef typename CFT::Numerator_type NUM;
            typedef typename CFT::Denominator_type DEN;

            if(ext.is_extended()){
                NUM a0_num, a1_num;
                DEN a0_den, a1_den;
                DEN common_den;
                decompose(ext.a0(),a0_num,a0_den);
                decompose(ext.a1(),a1_num,a1_den);
                common_den=common_factor(a0_den,a1_den);
                typename CGAL::Coercion_traits<NUM,DEN>::Cast cast;
                a0_num = cast(a0_num) * 
                         cast(CGAL::integral_division(a1_den,common_den));
                a1_num = cast(a1_num) * 
                         cast(CGAL::integral_division(a0_den,common_den)); 
                den = CGAL::integral_division(a0_den,common_den)*a1_den;
                num = Numerator_type(a0_num,a1_num,ext.root());
            }else{
                NUM a0_num;
                decompose(ext.a0(),a0_num,den);
                num = Numerator_type(a0_num);
            }
        }
    };
    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef NT result_type;
        NT operator () (const Numerator_type&   num,
                        const Denominator_type& den){
            if(num.is_extended()){
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                COEFF a1=compose(num.a1(),den);
                return NT(a0,a1,num.root());
            }else{
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                return NT(a0);
            }
        }
    };
};

template <class EXT, class BOOL_TAG>
class Sqrt_ext_Ftr_base_1;

template <class COEFF, class ROOT, class ACDE_TAG, class FP_TAG>
class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>, CGAL::Tag_true >
    : public Sqrt_ext_Ftr_base_2<
    Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>,
    ::boost::is_same< typename CGAL::Coercion_traits<ROOT,typename CGAL::Fraction_traits<COEFF>::Numerator_type>::Type,
                        typename CGAL::Fraction_traits<COEFF>::Numerator_type>::value >
{
    //nothing new
};

    template <class COEFF, class ROOT, class ACDE_TAG, class FP_TAG>
    class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>, CGAL::Tag_false >
      : public Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG>, false>
    {
        //nothing new
    };
} // namespace Intern


/*
namespace Intern{
    template <class SqrtExt,class BoolTag> class Sqrt_ext_Coftr_base_1;
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_false >{
    public:
        typedef SqrtExt          Numerator_type;
        typedef ::CGAL::Tag_false Is_composable;
        typedef ::CGAL::Null_tag Denominator_type;
        typedef ::CGAL::Null_tag Type;
        typedef ::CGAL::Null_tag Compose;
    };
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_true >{
        typedef typename SqrtExt::NT Coeff;
        typedef typename SqrtExt::ROOT Root;
        typedef typename CGAL::Cofraction_traits<Coeff> CFT;
        typedef typename CFT::Type Type_coeff;

    public:
        typedef SqrtExt                                       Numerator_type;
        typedef ::CGAL::Tag_true                               Is_composable;
        typedef typename CFT::Denominator_type                Denominator;
        typedef CGAL::Sqrt_extension<Type_coeff,Root,ACDE_TAG,FP_TAG> Type;

        class Compose {
    public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type
            typedef Denominator_type second_argument_type;
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type   num,
                                      Denominator_type den){
                if(num.is_extended()){
                    typename CFT::Compose compose_coeff;
                    Type_coeff a0_new(compose_coeff(num.a0(),den));
                    Type_coeff a1_new(compose_coeff(num.a1(),den));
                    return result_type(a0_new, a1_new, num.root());
                }else{
                    typename CFT::Compose compose_coeff;
                    return result_type(compose_coeff(num.a0(),den));
                }
            };
        };
    };
}

template <class Coeff, class Root,class ACDE_TAG, class FP_TAG>
class Cofraction_traits<Sqrt_extension<Coeff,Root,ACDE_TAG,FP_TAG> >
    :public Intern::Sqrt_ext_Coftr_base_1<
    Sqrt_extension<Coeff,Root,ACDE_TAG,FP_TAG>,
    typename CGAL::Cofraction_traits<Coeff>::Is_composable>{
    //nothing new;
};
*/

} //namespace CGAL

#endif
