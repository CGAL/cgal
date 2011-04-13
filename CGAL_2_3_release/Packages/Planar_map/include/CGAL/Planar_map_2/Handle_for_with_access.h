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

template <class RefCounted,
          class Allocator = CGAL_ALLOCATOR(RefCounted) >
class Handle_for_with_access : public Handle_for<RefCounted,Allocator>
{
  public:
  Handle_for_with_access(const RefCounted& rc) : 
    Handle_for<RefCounted,Allocator>(rc){}
  Handle_for_with_access() : Handle_for<RefCounted,Allocator>(){}
  Handle_for_with_access( const Handle_for_with_access& h) : 
    Handle_for<RefCounted,Allocator>(h){}
  
  typename Allocator::pointer pointer() const {return ptr;}
};

} // namespace CGAL
#endif // CGAL_PLANAR_MAP_2_HANDLE_FOR_WITH_ACCESS_H












