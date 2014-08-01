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

#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP
#include <iostream>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Gmpfr.h>
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
        typedef CGAL::Gmpfr     Gmpfr;
        typedef NumberType      NT;
        Gmpfr a(5.0);
        NT b(2);
        if((-a)==-5&&(a+b)==7&&(a-b)==3&&(a*b)==10&&(a/b)==2.5)
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

template<>
int test_operators<CGAL::Gmpfr>(){
        typedef CGAL::Gmpfr     Gmpfr;
        Gmpfr a(5.0);
        Gmpfr b(2);
        if((-a)==-5&&(a+b)==7&&(a-b)==3&&(a*b)==10&&(a/b)==2.5&&(a%b)==1)
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

template<class NumberType>
int test_arithmetic(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef NumberType      NT;
        Gmpfr a(7.0);
        NT b(2);
        Gmpfr::set_default_rndmode(std::round_to_nearest);
        if(
              CGAL::Gmpfr::add(a,b) == 9
           && CGAL::Gmpfr::add(a,b,std::round_to_nearest) == 9
           && CGAL::Gmpfr::add(a,b,2) == 8
           && CGAL::Gmpfr::add(a,b,2,std::round_toward_infinity) == 12
           && CGAL::Gmpfr::sub(a,b) == 5
           && CGAL::Gmpfr::sub(a,b,std::round_to_nearest) == 5
           && CGAL::Gmpfr::sub(a,b,2) == 4
           && CGAL::Gmpfr::sub(a,b,2,std::round_toward_infinity) == 6
           && CGAL::Gmpfr::mul(a,b) == 14
           && CGAL::Gmpfr::mul(a,b,std::round_to_nearest) == 14
           && CGAL::Gmpfr::mul(a,b,2) == 16
           && CGAL::Gmpfr::mul(a,b,2,std::round_toward_neg_infinity) == 12
           && CGAL::Gmpfr::div(a,b) == 3.5
           && CGAL::Gmpfr::div(a,b,std::round_to_nearest) == 3.5
           && CGAL::Gmpfr::div(a,b,2) == 4
           && CGAL::Gmpfr::div(a,b,2,std::round_toward_neg_infinity) == 3
           && Gmpfr(9.0).is_square(a) && a==3
        )
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

int test_to_integer_exp(CGAL::Gmpfr f){
        CGAL::Gmpq q;
        std::pair<CGAL::Gmpz,long> ze=f.to_integer_exp();
        if(ze.second>0)
                q=ze.first<<ze.second;
        else{
                q=ze.first;
                q/=CGAL::Gmpz(1)<<-(ze.second);
        }
        if(f==q)
                return 0;
        else{
                std::cerr<<"failed!"<<std::endl;
                exit(-1);
        }
}

// This function checks equality between an NumberType x and a Gmpfr y.
template<class NumberType>
int are_different(const NumberType &x,const CGAL::Gmpfr &y){
        return x!=y;
}

template<>
int are_different(const std::pair<CGAL::Gmpz,long> &x,const CGAL::Gmpfr &y){
        return(mpfr_cmp_si_2exp(y.fr(),mpz_get_si(x.first.mpz()),x.second));
}

template<class NumberType>
int test_constructors(const NumberType &x){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef NumberType      NT;
        CGAL_USE_TYPE(NT);
        bool fail=false;
        Gmpfr::set_default_precision(70);
        Gmpfr f(x);
        // this conversion should be exact
        if(are_different(x,f)){
                std::cerr<<"failed default construction! (inexact)"<<std::endl;
                fail=true;
        }
        Gmpfr g(x,100);
        if(g.get_precision()!=100){
                std::cerr<<"failed construction prec!"<<std::endl;
                fail=true;
        }
        Gmpfr h(x,std::round_toward_infinity);
        if(h.get_precision()!=Gmpfr::get_default_precision()){
                std::cerr<<"failed construction rndmode!"<<std::endl;
                fail=true;
        }
        Gmpfr i(x,std::round_toward_infinity,100);
        if(i.get_precision()!=100){
                std::cerr<<"failed construction prec/rndmode!"<<std::endl;
                fail=true;
        }
        if(fail)
                exit(-1);
        return 0;
}

int main(){
  typedef CGAL::Gmpfr NT;
  typedef CGAL::Field_with_kth_root_tag Tag;
  typedef CGAL::Tag_false Is_exact;

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

  TEST("constructors int",test_constructors<int>(-3);)
  TEST("constructors long",test_constructors<long>(-456);)
  TEST("constructors unsigned long",test_constructors<unsigned long>(78);)
  TEST("constructors double",test_constructors<double>(7.5);)
  TEST("constructors long double",test_constructors<long double>(-7.5);)
  TEST("constructors Gmpz",
       test_constructors<CGAL::Gmpz>((CGAL::Gmpz(1)<<1000)+CGAL::Gmpz(1));)
  TEST("constructors Gmpzf",test_constructors<CGAL::Gmpzf>(1025);)
  typedef std::pair<CGAL::Gmpz,long>                            MantExp;
  TEST("constructors pair<Gmpz,long>",
       test_constructors<MantExp>(std::make_pair(CGAL::Gmpz(4096),35));)

  TEST("operators Gmpfr",test_operators<NT>();)
  TEST("operators Gmpzf",test_operators<CGAL::Gmpzf>();)
  TEST("operators Gmpz",test_operators<CGAL::Gmpz>();)
  TEST("operators int",test_operators<int>();)
  TEST("operators long",test_operators<long>();)
  TEST("operators unsigned long",test_operators<unsigned long>();)

  TEST("arithmetic Gmpfr",test_arithmetic<NT>();)
  TEST("arithmetic long",test_arithmetic<long>();)
  TEST("arithmetic unsigned long",test_arithmetic<unsigned long>();)
  TEST("arithmetic int",test_arithmetic<int>();)
  TEST("arithmetic Gmpz",test_arithmetic<CGAL::Gmpz>();)
  TEST("arithmetic Gmpzf",test_arithmetic<CGAL::Gmpzf>();)

  NT plus_infinity;
  NT minus_infinity;
  CGAL::Gmpzf f;
  mpfr_set_inf(plus_infinity.fr(),1);
  mpfr_set_inf(minus_infinity.fr(),-1);
  TEST("I/O positive int",(CGAL::test_io<NT,int>(2145338339));)
  TEST("I/O negative int",(CGAL::test_io<NT,int>(-25295236));)
  TEST("I/O double",(CGAL::test_io<NT,double>(.2147483647));)
  TEST("I/O +inf",(CGAL::test_io<NT,NT>(plus_infinity));)
  TEST("I/O -inf",(CGAL::test_io<NT,NT>(minus_infinity));)
  TEST("I/O Gmpzf",(CGAL::test_io<NT,CGAL::Gmpzf>(f));)

  std::stringstream ss;
  CGAL::Gmpz z;
  ss<<"-4503599627370496";
  ss>>z;
  TEST("to_integer_exp e>0",(test_to_integer_exp(NT(std::make_pair(z,-50))));)
  z=CGAL::Gmpz(2147483647)*CGAL::Gmpz(2145338339);
  NT::Precision_type zsize=
          static_cast<NT::Precision_type>( mpz_sizeinbase(z.mpz(),2) );
  TEST("to_integer_exp e==0",(test_to_integer_exp(NT(z,zsize)));)
  TEST("to_integer_exp e<0",(test_to_integer_exp(NT(std::make_pair(z,-97))));)

  // TODO: missing tests for conversion functions
  // to_double, to_interval, to_double_exp, to_interval_exp

  return 0;
}

#undef TEST
#undef TEST_AS

#else
int main() { return 0; }
#endif
