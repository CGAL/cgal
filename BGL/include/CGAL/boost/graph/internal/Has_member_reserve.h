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
//
// Author(s) : Sebastien Loriot


#ifndef CGAL_HAS_MEMBER_RESERVE_H
#define CGAL_HAS_MEMBER_RESERVE_H

namespace CGAL {
namespace internal {

template<class T, class I1, class I2, class I3>
class Has_member_reserve
{
private:
  template<class U, U>
  class check {};

  template<class C>
  static char f(check<void(C::*)(I1, I2, I3), &C::reserve>*);

  template<class C>
  static int f(...);
public:
  static const bool value = (sizeof(f<T>(0)) == sizeof(char));
};

}  // internal
}  // CGAL

#endif /* CGAL_HAS_MEMBER_RESERVE_H */
