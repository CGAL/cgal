// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Philipp Moeller


#ifndef CGAL_HAS_MEMBER_CLEAR_H
#define CGAL_HAS_MEMBER_CLEAR_H

namespace CGAL {
namespace internal {

template<class T>
class Has_member_clear
{
private:
  template <class C>
  static auto f(int) -> decltype(std::declval<C>().clear(), char());

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

template<class T>
inline constexpr bool Has_member_clear_v = Has_member_clear<T>::value;

}  // internal
}  // cgal

#endif /* CGAL_HAS_MEMBER_CLEAR_H */
