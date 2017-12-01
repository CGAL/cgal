// Copyright (c) 2017 GeometryFactory (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_HAS_CONVERSION_H
#define CGAL_HAS_CONVERSION_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Homogeneous_converter.h>
#include <CGAL/representation_tags.h>

namespace CGAL {

namespace internal {

template<typename K1, typename K2, typename Rep = typename K1::Rep_tag /* Cartesian_tag */>
struct Converter_selector
{
  CGAL_static_assertion_msg((boost::is_same<typename K1::Rep_tag,
                                            typename K2::Rep_tag>::value),
                            "Kernels must have the same representation");

  typedef CGAL::Cartesian_converter<K1, K2> type;
};

template<typename K1, typename K2>
struct Converter_selector<K1, K2, Homogeneous_tag>
{
  CGAL_static_assertion_msg((boost::is_same<typename K1::Rep_tag,
                                            typename K2::Rep_tag>::value),
                            "Kernels must have the same representation");

  typedef CGAL::Homogeneous_converter<K1, K2> type;
};

} // namespace internal

/// Check whether Cartesian/Homogeneous_converter<K1, K2> provides either
/// of the following conversions:
///
///  -- K2T operator()(const K1T& ) const
///  -- const K2T& operator()(const K1T& ) const
///
/// \pre K1 and K2 are kernels with the same representation tag
/// \pre K1T is a typedef in K1
/// \pre K2T is a typedef in K2
template<typename K1, typename K2, typename K1T, typename K2T>
class Has_conversion
{
  typedef char one;
  typedef struct { char arr[2]; } two;

  template<typename U, U>
  struct Wrapper { };

  template<typename CC, typename T1, typename T2>
  static one test(Wrapper<T2 (CC::*)(const T1&) const, &CC::operator()>*);
  template<typename CC, typename T1, typename T2>
  static one test(Wrapper<const T2& (CC::*)(const T1&) const, &CC::operator()>*);

  template<typename CC, typename T1, typename T2>
  static two test(...);

public:
  static const bool value =
    sizeof(test<typename internal::Converter_selector<K1, K2>::type, K1T, K2T >(0)) == 1;
};

} // namespace CGAL

#endif // CGAL_HAS_CONVERSION_H
