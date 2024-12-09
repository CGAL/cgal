// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Jane Tournois


#ifndef CGAL_HAS_MEMBER_ID_H
#define CGAL_HAS_MEMBER_ID_H

namespace CGAL {
namespace internal {

  template <typename Type>
  class Has_member_id
  {
    template <typename U>
    static auto check(int) -> decltype(std::declval<U>().id(), char());

    template <typename U>
    static int check(...);

  public:
    static const bool value = (sizeof(char) == sizeof(check<Type>(0)));
  };

}  // internal
}  // cgal

#endif /* CGAL_HAS_MEMBER_ID_H */
