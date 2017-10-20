// Copyright (c) 2016 GeometryFactory (France).  All rights reserved.
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
// Author(s) : Jane Tournois


#ifndef CGAL_HAS_MEMBER_ID_H
#define CGAL_HAS_MEMBER_ID_H

namespace CGAL {
namespace internal {

  template <typename Type>
  class Has_member_id
  {
    typedef char yes[1];
    typedef char no[2];

    struct BaseWithId
    {
      void id(){}
    };
    struct Base : public Type, public BaseWithId {};

    template <typename T, T t>
    class Helper{};

    template <typename U>
    static no &check(U*, Helper<void (BaseWithId::*)(), &U::id>* = 0);

    static yes &check(...);

  public:
    static const bool value = (sizeof(yes) == sizeof(check((Base*)(0))));
  };

}  // internal
}  // cgal

#endif /* CGAL_HAS_MEMBER_ID_H */
