// Copyright (c) 1999  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Oren Nechushtan
 

#ifndef CGAL_PLANAR_MAP_2_HANDLE_FOR_WITH_ACCESS_H
#define CGAL_PLANAR_MAP_2_HANDLE_FOR_WITH_ACCESS_H

#ifndef CGAL_HANDLE_FOR_H
#include <CGAL/Handle_for.h>
#endif

namespace CGAL {

template <class T,
          class Allocator_ = CGAL_ALLOCATOR(T) >
class Handle_for_with_access : public Handle_for<T, Allocator_>
{
  public:
  Handle_for_with_access(const T& rc) : 
    Handle_for<T, Allocator_>(rc){}
  Handle_for_with_access() : Handle_for<T, Allocator_>(){}
  Handle_for_with_access( const Handle_for_with_access& h) : 
    Handle_for<T, Allocator_>(h){}
  
  const T* pointer() const {return Ptr();}

  T* pointer() {return ptr();}
  
};

} // namespace CGAL
#endif // CGAL_PLANAR_MAP_2_HANDLE_FOR_WITH_ACCESS_H












