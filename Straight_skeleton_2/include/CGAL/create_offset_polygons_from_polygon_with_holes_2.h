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
#ifndef CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H
#define CGAL_CREATE_OFFSET_POLYGONS_FROM_POLYGON_WITH_HOLES_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/disable_warnings.h>

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

#include <CGAL/enable_warnings.h>

#endif 
// EOF //
