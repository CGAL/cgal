// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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

// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/arrange_offset_polygons_2.h>

namespace CGAL {

template<class FT, class K>
std::vector< boost::shared_ptr< Polygon_with_holes_2<K> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2 ( FT const& aOffset, Polygon_with_holes_2<K> const& aPolyWH, K const& k = K() )
{
  return arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2( aOffset
                                                                                 , aPolyWH.outer_boundary()
                                                                                 , aPolyWH.holes_begin   ()
                                                                                 , aPolyWH.holes_end     ()
                                                                                 , k
                                                                                 )
                                 );
}

template<class FT, class K>
std::vector< boost::shared_ptr< Polygon_with_holes_2<K> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2 ( FT const& aOffset, Polygon_with_holes_2<K> const& aPolyWH, K const& k = K() )
{
  return arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_2( aOffset
                                                                                 , aPolyWH.outer_boundary()
                                                                                 , aPolyWH.holes_begin   ()
                                                                                 , aPolyWH.holes_end     ()
                                                                                 , k
                                                                                 )
                                 );
}


} //namespace CGAL


#endif 
// EOF //
