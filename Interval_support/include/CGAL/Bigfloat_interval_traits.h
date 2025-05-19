// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//                 Ron Wein         <wein@post.tau.ac.il>


#ifndef CGAL_BIGFLOAT_INTERVAL_TRAITS_H
#define CGAL_BIGFLOAT_INTERVAL_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>

namespace CGAL {

// TODO: rename this into MPFI_traits ?
// add a better rounding control

template<typename BigfloatInterval> class Bigfloat_interval_traits;

template<typename BFI> inline long get_significant_bits(BFI bfi) {
  typename Bigfloat_interval_traits<BFI>::Relative_precision relative_precision;
  return  zero_in(bfi) ? -1 : (std::max)(long(0),relative_precision(bfi));
}

template<typename BFI> inline long set_precision(BFI,long prec) {
    typename Bigfloat_interval_traits<BFI>::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BFI> inline long get_precision(BFI) {
    typename Bigfloat_interval_traits<BFI>::Get_precision get_precision;
    return get_precision();
}

template<typename BFI> inline long relative_precision(const BFI& bfi) {
    typename Bigfloat_interval_traits<BFI>::Relative_precision
      relative_precision;
    return relative_precision(bfi);
}
} //namespace CGAL

#endif // CGAL_BIGFLOAT_INTERVAL_TRAITS_H
