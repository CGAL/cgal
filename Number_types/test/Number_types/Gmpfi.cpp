// Copyright (c) 2009,2010 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifdef CGAL_USE_MPFI
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/_test_io.h>

#define TEST(_string,_code) \
        std::cerr<<"testing "<<_string<<": "<<std::flush; \
        _code; \
        std::cerr<<"OK"<<std::endl;

#define TEST_AS(_a,_b,_c) \
        std::cerr<<"testing algebraic structure ("<< \
                _a<<','<<_b<<','<<_c<<"): "<<std::flush; \
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(_a),NT(_b),NT(_c)); \
        std::cerr<<"OK"<<std::endl;

template<class NumberType>
int test_operators(){
        typedef CGAL::Gmpfi     Gmpfi;
        typedef NumberType      NT;
        Gmpfi a(5.0);
        NT b(2);
        if((-a)==-5&&(a+b)==7&&(a-b)==3&&(a*b)==10&&(a/b)==2.5)
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

template<class NumberType>
int test_arithmetic(){
        typedef CGAL::Gmpfi     Gmpfi;
        typedef NumberType      NT;
        Gmpfi a(7.0);
        NT b(2);
        if(
              CGAL::Gmpfi::add(a,b)            == 9
           && CGAL::Gmpfi::add(b,a,2).inf()    == 8
           && CGAL::Gmpfi::sub(a,b)            == 5
           && CGAL::Gmpfi::sub(b,a,2).sup()    == -4
           && CGAL::Gmpfi::mul(a,b)            == 14
           && CGAL::Gmpfi::mul(b,a,2).sup()    == 16
           && CGAL::Gmpfi::div(a,b)            == 3.5
           && CGAL::Gmpfi::div(b,a,2).inf()    == .25
           && Gmpfi(9.0).is_square(a) && a==3
        )
                return 0;
        else{
                std::cerr<<"error: arithmetic test"<<std::endl;
                exit(-1);
        }
}

int test_precision(){
        typedef CGAL::Gmpfi     Gmpfi;

        Gmpfi::Precision_type new_precision=200;
        Gmpfi::Precision_type old_precision;
        Gmpfi one(1);
        old_precision=Gmpfi::set_default_precision(new_precision);
        // this does not make sense when new_precision==old_precision
        Gmpfi two(2);
        if(old_precision==new_precision ||
           one.get_precision()!=old_precision ||
           two.get_precision()!=new_precision){
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
        return 0;
}

int test_refcount(){
        CGAL::Gmpq A("6420587669/17179869184");
        CGAL::Gmpfi f(A,20);
        CGAL::Gmpfr i=f.inf();
        assert(!i.is_unique());
        i*=2;
        assert(i!=f.inf());
        assert(i.is_unique());
        return 0;
}

int main(){
        typedef CGAL::Gmpfi                     NT;
        typedef CGAL::Field_with_kth_root_tag   Tag;
        typedef CGAL::Tag_false                 Is_exact;

        TEST("algebraic structure (default)",
             (CGAL::test_algebraic_structure<NT,Tag,Is_exact>()))
        TEST_AS(4,6,15)
        TEST_AS(-4,6,15)
        TEST_AS(4,-6,15)
        TEST_AS(-4,-6,15)
        TEST_AS(4,6,-15)
        TEST_AS(-4,6,-15)
        TEST_AS(4,-6,-15)
        TEST_AS(-4,-6,-15)
        TEST("real embeddable",CGAL::test_real_embeddable<NT>())

        TEST("operators Gmpfi",test_operators<NT>())
        TEST("operators Gmpfr",test_operators<CGAL::Gmpfr>();)
        TEST("operators Gmpz",test_operators<CGAL::Gmpz>();)
        TEST("operators Gmpq",test_operators<CGAL::Gmpq>();)
        TEST("operators int",test_operators<int>();)
        TEST("operators long",test_operators<long>();)
        TEST("operators unsigned long",test_operators<unsigned long>();)

        TEST("arithmetic Gmpfi",test_arithmetic<CGAL::Gmpfi>();)
        TEST("arithmetic Gmpfr",test_arithmetic<CGAL::Gmpfr>();)
        TEST("arithmetic long",test_arithmetic<long>();)
        TEST("arithmetic unsigned long",test_arithmetic<unsigned long>();)
        TEST("arithmetic int",test_arithmetic<int>();)
        TEST("arithmetic Gmpz",test_arithmetic<CGAL::Gmpz>();)
        TEST("arithmetic Gmpq",test_arithmetic<CGAL::Gmpq>();)

        CGAL::Gmpfr plus_infinity;
        CGAL::Gmpfr minus_infinity;
        CGAL::Gmpzf zf;

        mpfr_set_inf(plus_infinity.fr(),1);
        mpfr_set_inf(minus_infinity.fr(),-1);

        TEST("I/O int",(CGAL::test_interval_io<NT,int>(2145338339)))
        TEST("I/O double",(CGAL::test_interval_io<NT,double>(.2147483647)))
        TEST("I/O Gmpfr",
             (CGAL::test_interval_io<NT,CGAL::Gmpfr>
              (CGAL::Gmpfr(1,100)/CGAL::Gmpfr(3,100))))
        TEST("I/O [int,+inf]",
             (CGAL::test_interval_io<NT,std::pair<int,CGAL::Gmpfr> >
              (std::make_pair(2145338339,plus_infinity))))
        TEST("I/O [-inf,double]",
             (CGAL::test_interval_io<NT,std::pair<CGAL::Gmpfr,double> >
              (std::make_pair(minus_infinity,.2147483647))))
        TEST("I/O [-inf,+inf]",
             (CGAL::test_interval_io<NT,std::pair<CGAL::Gmpfr,CGAL::Gmpfr> >
              (std::make_pair(minus_infinity,plus_infinity))))

        TEST("precision",test_precision())

        // TODO: missing tests for conversion functions
        // to_double, to_interval, to_double_exp, to_interval_exp

#ifndef CGAL_GMPFR_NO_REFCOUNT
        TEST("endpoint reference counting",test_refcount())
#else
        std::cerr<<"endpoint reference counting was not tested"<<std::endl;
#endif

        return 0;
}

#undef TEST
#undef TEST_AS

#else
int main(){
        return 0;
}
#endif
