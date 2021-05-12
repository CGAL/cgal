// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_DEFAULT_H
#define CGAL_HALFEDGEDS_DEFAULT_H 1

#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/memory.h>
#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>


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
