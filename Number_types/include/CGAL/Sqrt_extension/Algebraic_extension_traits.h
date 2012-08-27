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


#ifndef CGAL_SQRT_EXTENSION_ALGEBRAIC_EXTENSION_TRAITS_H
#define CGAL_SQRT_EXTENSION_ALGEBRAIC_EXTENSION_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
class Algebraic_extension_traits<CGAL::Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> > {
/* needed to 'add up' sqrt_extensions in iterator range such that all roots 
   are collected in order to keep operation time minimal all scalar coeffs 
   are set to 1 by standardise. 
   TODO .. find a better name, if you want to.
*/
    template <class Type_>
    class Standardise {
    public:
        typedef Type_ argument_type;
        typedef Type_ result_type;
    Type_ operator () (const Type_&) const {
            return Type_(1);
        }
    };
    
  template <class COEFF_, class ROOT_,class ACDE_TAG_,class FP_TAG_>
  class Standardise<CGAL::Sqrt_extension<COEFF_,ROOT_,ACDE_TAG_,FP_TAG_> > {
        Standardise<COEFF_> standardise;
    public:
    typedef CGAL::Sqrt_extension<COEFF_,ROOT_,ACDE_TAG_,FP_TAG_> Type_;
        typedef Type_ argument_type;
        typedef Type_ result_type;
    Type_ operator () (const Type_& a) const {       
            if (a.a1() != COEFF_(0)){
                return Type_(standardise(a.a0()),standardise(a.a1()),a.root());
            }else{
                return Type_(standardise(a.a0()));
            }
        }
    };

public:
    //! \name Typedefs 
    //! the number type for which this instance has been instantiated
    typedef CGAL::Sqrt_extension<COEFF,ROOT,ACDE_TAG,FP_TAG> Type;
    //! Sqrt_extension as a number type is extended
    typedef ::CGAL::Tag_true Is_extended;
    
    /*! computes the factor which normalizes a number to be integral 
     *  after multiplication
     */
    class Normalization_factor{
    public:
        //! argument type
        typedef Type argument_type;
        //! result type
        typedef Type result_type;
        //! determine normalization factor
        Type operator () (const Type& a) {
            typedef Algebraic_extension_traits<COEFF> SET;
            typename  SET::Normalization_factor normalization_factor;
            CGAL_precondition(a != Type(0));

            Type result;
            if(a.is_extended() && a.a1() != COEFF(0) ){
                Type tmp1(a.a0(),-a.a1(),a.root());
                Type tmp2= a*tmp1;
                CGAL_postcondition(tmp2.a1()==COEFF(0));
                result = tmp1*Type(normalization_factor(tmp2.a0()));
                CGAL_postcondition(result.a1() != COEFF(0));
            }else{
                result = Type(normalization_factor(a.a0())); 
            }           
            return result;
        }
    };

    /*! returns the extension factor needed for the gcd_utcf computation 
     *  for more details see ... TODO!!
     */
    class Denominator_for_algebraic_integers {
    public:
        //! argument type
        typedef Type argument_type;
        //! result type
        typedef Type result_type;
        //! determine denominator for algebraic integers
    public:
        Type operator () (const Type& a) {
            
            typedef Algebraic_extension_traits<COEFF> AET;
            typename AET::Denominator_for_algebraic_integers dfai;

            Standardise<COEFF> standardise;
            if (a.a1() != COEFF(0)) {
                COEFF tmp = 
                    standardise(a.a0())
                    + standardise(a.a1())
                    + standardise(COEFF(a.root()));
                return  Type(COEFF(4) * COEFF(a.root()))* Type(dfai(tmp));
            } else {
                return Type(dfai(a.a0()));
            };
        }

        /*! overloaded operator for computing the denominator out of an iterator
          range accumulates all root expressions which appear in the range to 
          the most complex term and uses this term to determine the denominator 
          for algebraic integers
        */
        template <class InputIterator>
        Type operator () (InputIterator begin, InputIterator end) {
            Type a = std::accumulate(::boost::make_transform_iterator(begin,Standardise<Type>()), 
                    ::boost::make_transform_iterator(end  ,Standardise<Type>()), Type(0));
            return (*this)(a);
        }
    };
};

} //namespace CGAL

#endif
