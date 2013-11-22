// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

// This file contains the simplest refiner, that bisects the interval a given
// number of times.

#ifndef CGAL_RS_BISECTION_REFINER_1_H
#define CGAL_RS_BISECTION_REFINER_1_H

#include <CGAL/Polynomial_traits_d.h>
#include "signat_1.h"

namespace CGAL{

template <class Polynomial_,class Bound_>
struct Bisection_refiner_1{
        typedef CGAL::RS_AK1::Signat_1<Polynomial_,Bound_>           Signat;
        void operator()(const Polynomial_&,Bound_&,Bound_&,int);
}; // class Bisection_refiner_1

// TODO: Write in a generic way, if possible.
template <class Polynomial_,class Bound_>
void
Bisection_refiner_1<Polynomial_,Bound_>::
operator()(const Polynomial_&,Bound_&,Bound_&,int){
        CGAL_error_msg("bisection refiner not implemented for these types");
        return;
}

// TODO: Optimize this function.
template<>
void
Bisection_refiner_1<Polynomial<Gmpz>,Gmpfr>::
operator()(const Polynomial<Gmpz> &pol,Gmpfr &left,Gmpfr &right,int prec){
        typedef Polynomial<Gmpz>                        Polynomial;
        typedef Polynomial_traits_d<Polynomial>         Ptraits;
        typedef Ptraits::Make_square_free               Sfpart;
        typedef CGAL::RS_AK1::Signat_1<Polynomial,Gmpfr>
                                                        Signat;
        CGAL_precondition(left<=right);
        // TODO: add precondition to check whether the interval is a point
        // or the evaluations on its endpoints have different signs
        //std::cout<<"refining ["<<left<<","<<right<<"]"<<std::endl;

#ifndef CGAL_GMPFR_NO_REFCOUNT
        // Make sure the endpoints do not share references. If some of them
        // does, copy it.
        if(!left.is_unique()){
                Gmpfr new_left(0,left.get_precision());
                mpfr_set(new_left.fr(),left.fr(),GMP_RNDN);
                left=new_left;
                CGAL_assertion_code(new_left=Gmpfr();)
                CGAL_assertion(left.is_unique());
        }
        if(!right.is_unique()){
                Gmpfr new_right(0,right.get_precision());
                mpfr_set(new_right.fr(),right.fr(),GMP_RNDN);
                right=new_right;
                CGAL_assertion_code(new_right=Gmpfr();)
                CGAL_assertion(right.is_unique());
        }
#endif // CGAL_GMPFR_NO_REFCOUNT

        Polynomial sfpp=Sfpart()(pol);
        Signat signof(sfpp);
        CGAL::Sign sl,sc;
        mp_prec_t pl,pc;
        mpfr_t center;

        sl=signof(left);
        if(sl==ZERO)
                return;
        pl=left.get_precision();
        pc=right.get_precision();
        pc=(pl>pc?pl:pc)+(mp_prec_t)prec;
        mpfr_init2(center,pc);
        CGAL_assertion_code(int round=)
        mpfr_prec_round(left.fr(),pc,GMP_RNDN);
        CGAL_assertion(!round);
        CGAL_assertion_code(round=)
        mpfr_prec_round(right.fr(),pc,GMP_RNDN);
        CGAL_assertion(!round);
        for(int i=0;i<prec;++i){
                CGAL_assertion_code(round=)
                mpfr_add(center,left.fr(),right.fr(),GMP_RNDN);
                CGAL_assertion(!round);
                CGAL_assertion_code(round=)
                mpfr_div_2ui(center,center,1,GMP_RNDN);
                CGAL_assertion(!round);
                sc=signof(Gmpfr(center));
                if(sc==ZERO){   // we have a root
                        CGAL_assertion_code(round=)
                        mpfr_set(left.fr(),center,GMP_RNDN);
                        CGAL_assertion(!round);
                        mpfr_swap(right.fr(),center);
                        break;
                }
                if(sc==sl)
                        mpfr_swap(left.fr(),center);
                else
                        mpfr_swap(right.fr(),center);
        }
        mpfr_clear(center);
        CGAL_postcondition(left<=right);
        //std::cout<<"ref root is ["<<left<<","<<right<<"]"<<std::endl;
        return;
}

} // namespace CGAL

#endif // CGAL_RS_BISECTION_REFINER_1_H
