// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Michael Hemmer <mhemmer@uni-mainz.de>
// 

// Test program for the CGAL::Root_of_traits 



#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Root_of_traits.h>
#include <CGAL/Lazy_exact_nt.h>  
#include <CGAL/Test/test_root_of_traits.h>




template <class Integer, class Rational> 
void test_root_of_traits_for_set(Integer, Rational, CGAL::Null_tag){
  {
    typedef Integer RT;
    typedef Rational FT;
    typedef CGAL::Sqrt_extension<FT,FT,CGAL::Tag_true,CGAL::Tag_true> Root_of_2;
    
    CGAL::Test::test_root_of_traits<RT,FT,Root_of_2>();
    CGAL::Test::test_root_of_traits<FT,FT,Root_of_2>();
  }{
    typedef CGAL::Lazy_exact_nt<Integer>  RT;
    typedef CGAL::Lazy_exact_nt<Rational> FT;
    typedef CGAL::Lazy_exact_nt<CGAL::Sqrt_extension<Rational,Rational,CGAL::Tag_true,CGAL::Tag_true> > Root_of_2;
    CGAL::Test::test_root_of_traits<RT,FT,Root_of_2>(); 
    CGAL::Test::test_root_of_traits<FT,FT,Root_of_2>();
  }
}

template <class Integer, class Rational, class FWS> 
void test_root_of_traits_for_set(Integer, Rational, FWS){
  CGAL::Test::test_root_of_traits<FWS,FWS,FWS>();
  test_root_of_traits_for_set(Integer(), Rational(), CGAL::Null_tag());
}

int main(){
    CGAL::Test::test_root_of_traits< double , double , double >();
    
#ifdef CGAL_USE_GMP
    //TODO: switch to Gmpq
    {
      typedef CGAL::GMP_arithmetic_kernel AK; 
      typedef AK::Integer Integer;
      typedef AK::Rational Rational;
      typedef AK::Field_with_sqrt FWS;
      test_root_of_traits_for_set(Integer(),Rational(),FWS());
    }
#endif
#ifdef CGAL_USE_LEDA
    //TODO: switch to Gmpq
    {
      typedef CGAL::LEDA_arithmetic_kernel AK; 
      typedef AK::Integer Integer;
      typedef AK::Rational Rational;
      typedef AK::Field_with_sqrt FWS;
      test_root_of_traits_for_set(Integer(),Rational(),FWS());
    }
#endif
#ifdef CGAL_USE_CORE
    //TODO: switch to Gmpq
    {
      typedef CGAL::CORE_arithmetic_kernel AK; 
      typedef AK::Integer Integer;
      typedef AK::Rational Rational;
      typedef AK::Field_with_sqrt FWS;
      test_root_of_traits_for_set(Integer(),Rational(),FWS());
    }
#endif
    {
      typedef CGAL::MP_Float_arithmetic_kernel AK; 
      typedef AK::Integer Integer;
      typedef AK::Rational Rational;
      typedef AK::Field_with_sqrt FWS;
      test_root_of_traits_for_set(Integer(),Rational(),FWS());
    }    
    
    return 0;
}
