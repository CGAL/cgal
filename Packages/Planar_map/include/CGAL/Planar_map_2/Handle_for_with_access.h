// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : $CGAL_Revision: CGAL-2.2-I-26 $
// release_date  : $CGAL_Date: 2000/07/11 $
// 
// file          : include/CGAL/Planar_map_2/Handle_for_with_access.h
// package       : Planar_map
// maintainer    : Oren Nechushtan <theoren@math.tau.ac.il>
// revision      : 1.0
// revision_date : 27 Jun 2000 
// author(s)     : Oren Nechushtan
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
//
// ======================================================================
 

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












