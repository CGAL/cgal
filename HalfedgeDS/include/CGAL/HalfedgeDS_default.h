// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_DEFAULT_H
#define CGAL_HALFEDGEDS_DEFAULT_H 1

#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/memory.h>

namespace CGAL {

template <class Traits_, class HalfedgeDSItems = HalfedgeDS_items_2, 
          class Alloc = CGAL_ALLOCATOR(int)>
class HalfedgeDS_default 
    : public HalfedgeDS_list< Traits_, HalfedgeDSItems, Alloc> {
public:
    typedef Traits_                                          Traits;
    typedef HalfedgeDS_list<Traits_, HalfedgeDSItems, Alloc> D_S;
    typedef typename D_S::size_type                           size_type;
    HalfedgeDS_default() {}
    HalfedgeDS_default( size_type v, size_type h, size_type f)
        : HalfedgeDS_list< Traits_, HalfedgeDSItems, Alloc>(v,h,f) {}
};
#define CGAL_HALFEDGEDS_DEFAULT  ::CGAL::HalfedgeDS_default

} //namespace CGAL
#endif // CGAL_HALFEDGEDS_DEFAULT_H //
// EOF //
