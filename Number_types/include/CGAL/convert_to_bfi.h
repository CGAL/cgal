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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>



#ifndef CGAL_CONVERT_TO_BFI_H
#define CGAL_CONVERT_TO_BFI_H

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>

CGAL_BEGIN_NAMESPACE

template <class NTX>
typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval BFI; 
    typedef CGAL::Coercion_traits<NTX,BFI> CT;
    return typename CT::Cast()(x);
    
    // typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AT;
    // typedef typename AT::Bigfloat_interval BFI;
    // typename Bigfloat_interval_traits<BFI>::Convert_to_bfi convert_to_bfi;
    // return convert_to_bfi(x);
}


// TODO: move this to sqrt_extension ?
/*
  template <typename A,typename B> class Sqrt_extension;
  template <typename NT, typename ROOT>
  typename Get_arithmetic_kernel<NT>::Arithmetic_kernel::Bigfloat_interval
  convert_to_bfi(const CGAL::Sqrt_extension<NT,ROOT>& x) {
  typedef typename Get_arithmetic_kernel<NT>::Arithmetic_kernel AT;
  typedef typename AT::Bigfloat_interval BFI;
  if(x.is_extended()){
  BFI a0(convert_to_bfi(x.a0()));
  BFI a1(convert_to_bfi(x.a1()));
  BFI root(CGAL::sqrt(convert_to_bfi(x.root())));
  return a0+a1*root;
  }else{
  return convert_to_bfi(x.a0());
  }
  }
*/

CGAL_END_NAMESPACE 

#endif // CGAL_CONVERT_TO_BFI_H
