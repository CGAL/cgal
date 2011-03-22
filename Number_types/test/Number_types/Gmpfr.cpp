// Copyright (c) 2009,2010 Inria Lorraine (France). All rights reserved.
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

#define _TEST(_string,_code) \
        std::cerr<<"testing "<<_string<<": "<<std::flush; \
        _code; \
        std::cerr<<"OK"<<std::endl;

#define _TEST_AS(_a,_b,_c) \
        std::cerr<<"testing algebraic structure ("<< \
                _a<<','<<_b<<','<<_c<<"): "<<std::flush; \
        CGAL::test_algebraic_structure<NT,Tag,Is_exact>(NT(_a),NT(_b),NT(_c)); \
        std::cerr<<"OK"<<std::endl;

template<class _NT>
int test_operators(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef _NT             NT;
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

template<class _NT>
int test_arithmetic(){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef _NT             NT;
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

template<class _NT>
int test_constructors(_NT x){
        typedef CGAL::Gmpfr     Gmpfr;
        typedef _NT             NT;
        bool fail=false;
        Gmpfr::set_default_precision(70);
        Gmpfr f(x);
        // this conversion should be exact
        if(f!=x){
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

  _TEST("algebraic structure (default)",
        (CGAL::test_algebraic_structure<NT,Tag,Is_exact>()))
  _TEST_AS(4,6,15)
  _TEST_AS(-4,6,15)
  _TEST_AS(4,-6,15)
  _TEST_AS(-4,-6,15)
  _TEST_AS(4,6,-15)
  _TEST_AS(-4,6,-15)
  _TEST_AS(4,-6,-15)
  _TEST_AS(-4,-6,-15)
  _TEST("real embeddable",CGAL::test_real_embeddable<NT>())

  _TEST("constructors int",test_constructors<int>(-3);)
  _TEST("constructors long",test_constructors<long>(-456);)
  _TEST("constructors unsigned long",test_constructors<unsigned long>(78);)
  _TEST("constructors double",test_constructors<double>(7.5);)
  _TEST("constructors long double",test_constructors<long double>(-7.5);)
  _TEST("constructors Gmpz",
        test_constructors<CGAL::Gmpz>((CGAL::Gmpz(1)<<1000)+CGAL::Gmpz(1));)
  _TEST("constructors Gmpzf",test_constructors<CGAL::Gmpzf>(1025);)

  _TEST("operators Gmpfr",test_operators<NT>();)
  _TEST("operators Gmpzf",test_operators<CGAL::Gmpzf>();)
  _TEST("operators Gmpz",test_operators<CGAL::Gmpz>();)
  _TEST("operators int",test_operators<int>();)
  _TEST("operators long",test_operators<long>();)
  _TEST("operators unsigned long",test_operators<unsigned long>();)

  _TEST("arithmetic Gmpfr",test_arithmetic<NT>();)
  _TEST("arithmetic long",test_arithmetic<long>();)
  _TEST("arithmetic unsigned long",test_arithmetic<unsigned long>();)
  _TEST("arithmetic int",test_arithmetic<int>();)
  _TEST("arithmetic Gmpz",test_arithmetic<CGAL::Gmpz>();)
  _TEST("arithmetic Gmpzf",test_arithmetic<CGAL::Gmpzf>();)

  NT plus_infinity;
  NT minus_infinity;
  CGAL::Gmpzf f;
  mpfr_set_inf(plus_infinity.fr(),1);
  mpfr_set_inf(minus_infinity.fr(),-1);
  _TEST("I/O int",(CGAL::test_io<NT,int>(2145338339));)
  _TEST("I/O double",(CGAL::test_io<NT,double>(.2147483647));)
  _TEST("I/O +inf",(CGAL::test_io<NT,NT>(plus_infinity));)
  _TEST("I/O -inf",(CGAL::test_io<NT,NT>(minus_infinity));)
  _TEST("I/O Gmpzf",(CGAL::test_io<NT,CGAL::Gmpzf>(f));)

  std::stringstream ss;
  CGAL::Gmpz z;
  ss<<"-4503599627370496";
  ss>>z;
  _TEST("to_integer_exp e>0",(test_to_integer_exp(NT(std::make_pair(z,-50))));)
  z=CGAL::Gmpz(2147483647)*CGAL::Gmpz(2145338339);
  NT::Precision_type zsize=mpz_sizeinbase(z.mpz(),2);
  _TEST("to_integer_exp e==0",(test_to_integer_exp(NT(z,zsize)));)
  _TEST("to_integer_exp e<0",(test_to_integer_exp(NT(std::make_pair(z,-97))));)

  // TODO: missing tests for conversion functions
  // to_double, to_interval, to_double_exp, to_interval_exp

  return 0;
}

#undef _TEST
#undef _TEST_AS

#else
int main() { return 0; }
#endif
