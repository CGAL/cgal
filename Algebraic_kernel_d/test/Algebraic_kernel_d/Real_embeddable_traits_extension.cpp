// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL:$
// $Id:$
// 
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief test for Real_embeddable_traits_extension



#include <CGAL/basic.h> 
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>

//#include <CGAL/Sqrt_extension.h> // used in this file


void test_real_embeddable_extension(CGAL::Null_tag){ return ; }

template <typename NT_>
void test_real_embeddable_extension(const NT_&){
  typedef typename CGAL::Get_arithmetic_kernel<NT_>::Arithmetic_kernel AK;
  typedef typename AK::Integer Integer;

  typedef CGAL::internal::Real_embeddable_extension<NT_> RETE;

  typedef typename RETE::Type NT; 
  
  typedef typename RETE::Floor Floor;
  typedef typename RETE::Ceil Ceil; 
  typedef typename RETE::Floor_log2_abs Floor_log2_abs;
  typedef typename RETE::Ceil_log2_abs Ceil_log2_abs;
  
  {
    const Floor floor = Floor();
    typedef typename Floor::argument_type Argument_type;
    typedef typename Floor::result_type   Result_type;
    CGAL_static_assertion(( ::boost::is_same<NT, Argument_type>::value));  
    CGAL_static_assertion(( ::boost::is_same<Integer, Result_type>::value));
    assert(Integer(42) == floor(NT(42)));
    assert(Integer(-42) == floor(NT(-42)));
  }

  {
    const Floor_log2_abs floor_log2_abs = Floor_log2_abs();
    typedef typename Floor_log2_abs::argument_type Argument_type;
    typedef typename Floor_log2_abs::result_type   Result_type;
    CGAL_static_assertion(( ::boost::is_same<NT, Argument_type>::value));  
    CGAL_static_assertion(( ::boost::is_same<long, Result_type>::value));
    
    assert(long(0) == floor_log2_abs(NT(1)));
    assert(long(0) == floor_log2_abs(NT(-1)));

    assert(long(1) == floor_log2_abs(NT(2)));
    
    assert(long(1) == floor_log2_abs(NT(3)));
    assert(long(2) == floor_log2_abs(NT(4)));
    assert(long(2) == floor_log2_abs(NT(5)));
    
    assert(long(2) == floor_log2_abs(NT(7)));
    assert(long(3) == floor_log2_abs(NT(8)));
    assert(long(3) == floor_log2_abs(NT(9)));

    assert(long(2) == floor_log2_abs(NT(-7)));
    assert(long(3) == floor_log2_abs(NT(-8)));
    assert(long(3) == floor_log2_abs(NT(-9)));
  }

  {
    const Ceil ceil = Ceil();
    typedef typename Ceil::argument_type Argument_type;
    typedef typename Ceil::result_type   Result_type;
    CGAL_static_assertion(( ::boost::is_same<NT, Argument_type>::value));  
    CGAL_static_assertion(( ::boost::is_same<Integer, Result_type>::value));
    assert(Integer(42) == ceil(NT(42)));
    assert(Integer(-42) == ceil(NT(-42)));
  }

  {
    const Ceil_log2_abs ceil_log2_abs = Ceil_log2_abs();
    typedef typename Ceil_log2_abs::argument_type Argument_type;
    typedef typename Ceil_log2_abs::result_type   Result_type;
    CGAL_static_assertion(( ::boost::is_same<NT, Argument_type>::value));  
    CGAL_static_assertion(( ::boost::is_same<long, Result_type>::value));
    
    assert(long(0) == ceil_log2_abs(NT(1)));
    assert(long(0) == ceil_log2_abs(NT(-1)));

    assert(long(1) == ceil_log2_abs(NT(2)));
    
    assert(long(2) == ceil_log2_abs(NT(3)));
    assert(long(2) == ceil_log2_abs(NT(4)));
    assert(long(3) == ceil_log2_abs(NT(5)));
    
    assert(long(3) == ceil_log2_abs(NT(7)));
    assert(long(3) == ceil_log2_abs(NT(8)));
    assert(long(4) == ceil_log2_abs(NT(9)));

    assert(long(3) == ceil_log2_abs(NT(-7)));
    assert(long(3) == ceil_log2_abs(NT(-8)));
    assert(long(4) == ceil_log2_abs(NT(-9)));
  }
   
}



template <typename AK>
void test_real_embeddable_extension_ak(){
  typedef typename AK::Integer Integer;
  typedef typename AK::Bigfloat Bigfloat; 
  typedef typename AK::Bigfloat Bigfloat_interval; 

  test_real_embeddable_extension(Integer());
  //typedef typename AK::Rational Rational;
  //test_real_embeddable_extension(Rational()); TODO
  test_real_embeddable_extension(Bigfloat());
  test_real_embeddable_extension(Bigfloat_interval());
}


int main() {
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
    test_real_embeddable_extension_ak< CGAL::GMP_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
    test_real_embeddable_extension_ak< CGAL::LEDA_arithmetic_kernel >();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
    test_real_embeddable_extension_ak< CGAL::CORE_arithmetic_kernel >();
#endif
    
    return 0;
}
