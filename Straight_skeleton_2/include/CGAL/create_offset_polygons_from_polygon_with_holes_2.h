// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/arrange_offset_polygons_2.h>

namespace CGAL {

template<class FT, class OfK, class SsK, class C>
std::vector< boost::shared_ptr< Polygon_2<OfK,C> > >
inline
create_interior_skeleton_and_offset_polygons_2 ( FT const& aOffset, Polygon_with_holes_2<OfK,C> const& aPoly, SsK const& ssk )
{
  OfK ofk ;
  return create_interior_skeleton_and_offset_polygons_2(aOffset
                                                       ,aPoly.outer_boundary()
                                                       ,aPoly.holes_begin()
                                                       ,aPoly.holes_end()
                                                       ,ofk
                                                       ,ssk
                                                       );
    
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_2<OfK,C> > >
inline
create_interior_skeleton_and_offset_polygons_2 ( FT const& aOffset, Polygon_with_holes_2<OfK,C> const& aPoly )
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, Exact_predicates_inexact_constructions_kernel() );
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK,C> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2 ( FT const&             aOffset
                                                          , Polygon_2<OfK,C> const& aPoly
                                                          , bool                  aDontReverseOrientation = false
                                                          )
{
  return arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset
                                                                                            ,aPoly
                                                                                            ,aDontReverseOrientation
                                                                                            )
                                  );
}

  template<class FT, class OfK, class SsK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK,C> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2 ( FT const& aOffset, Polygon_with_holes_2<OfK,C> const& aPoly, SsK const& ssk )
{
  return arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(aOffset,aPoly,ssk));
}


  template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK,C> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2 ( FT const& aOffset, Polygon_with_holes_2<OfK,C> const& aPoly )
{
  return arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(aOffset,aPoly));
}

template<class FT, class OfK, class SsK>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2 ( FT const&             aOffset
                                                          , Polygon_2<OfK> const& aPoly
                                                          , SsK const&            ssk 
                                                          )
{
  return arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset
                                                                                            ,aPoly
                                                                                            ,ssk
                                                                                            )
                                  );
}



} // end namespace CGAL


#endif 
// EOF //
