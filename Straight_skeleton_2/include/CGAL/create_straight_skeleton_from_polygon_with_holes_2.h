// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//

// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CREATE_STRAIGHT_SKELETON_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_STRAIGHT_SKELETON_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

template<class K, class C>
boost::shared_ptr< Straight_skeleton_2<K> >
inline
create_interior_straight_skeleton_2 ( Polygon_with_holes_2<K,C> const& aPolyWithHoles )
{
  return create_interior_straight_skeleton_2(aPolyWithHoles.outer_boundary().vertices_begin()
                                            ,aPolyWithHoles.outer_boundary().vertices_end  ()
                                            ,aPolyWithHoles.holes_begin   ()
                                            ,aPolyWithHoles.holes_end     ()
                                            ,K()
                                            );
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_2_H //
// EOF //
