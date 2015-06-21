// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany),
// INRIA Sophia-Antipolis (France). 
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//
//    \brief provides tests convert_to_bfi.


#include <CGAL/basic.h>

#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

#include <cassert>
#include <CGAL/tags.h>
#include <CGAL/use.h>

#include <CGAL/convert_to_bfi.h> 
#include <CGAL/Sqrt_extension.h> 

#ifndef CGAL_TEST_CONVERT_TO_BFI_H
#define CGAL_TEST_CONVERT_TO_BFI_H

namespace CGAL {

// BFI = Bigfloat_interval 

template <typename BFI>
void test_convert_to_bfi_from(BFI,CGAL::Null_tag){return;}

template <typename BFI, typename From>
void test_convert_to_bfi_from(BFI,From){
  typedef typename CGAL::Coercion_traits<BFI,From>::Type CT_type;
  CGAL_USE_TYPE(CT_type);
  CGAL_static_assertion(( ::boost::is_same<CT_type, BFI>::value));
  assert(CGAL::convert_to_bfi(From(0))  == BFI(0));
  assert(CGAL::convert_to_bfi(From(1))  == BFI(1));
  assert(CGAL::convert_to_bfi(From(2))  == BFI(2));
  assert(CGAL::convert_to_bfi(From(4))  == BFI(4));
  assert(CGAL::convert_to_bfi(From(-8)) == BFI(-8));
  return;
}

template <typename BFI>
void test_convert_to_bfi(){
  typedef typename CGAL::Get_arithmetic_kernel<BFI>::Arithmetic_kernel AK;
  typedef typename AK::Integer Integer; 
  typedef typename AK::Rational Rational; 
  typedef typename AK::Field_with_sqrt Field_with_sqrt; // may be Null_tag 
  typedef CGAL::Sqrt_extension<Integer,Integer> EXT;

  test_convert_to_bfi_from( BFI() , Integer());
  test_convert_to_bfi_from( BFI() , Rational());
  test_convert_to_bfi_from( BFI() , Field_with_sqrt());
  test_convert_to_bfi_from( BFI() , EXT());
}

} //namespace CGAL

#endif // CGAL_TEST_REAL_COMPARABLE_H
