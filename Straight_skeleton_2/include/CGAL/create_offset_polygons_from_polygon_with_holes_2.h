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

#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/arrange_offset_polygons_2.h>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

// create_interior_skeleton_and_offset_polygons_2

template<class FT, class OfK, class SsK, class C>
std::vector< boost::shared_ptr< Polygon_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const Polygon_with_holes_2<OfK, C>& aPoly,
                                               const OfK& ofk,
                                               const SsK& ssk)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly.outer_boundary(),
                                                        aPoly.holes_begin(), aPoly.holes_end(),
                                                        ofk, ssk);
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const Polygon_with_holes_2<OfK, C>& aPoly,
                                               const OfK& ofk)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ofk);
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const Polygon_with_holes_2<OfK, C>& aPoly)
{
  OfK ofk;
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ofk);
}

// create_interior_skeleton_and_offset_polygons_with_holes_2

template<class FT, class OfK, class SsK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_with_holes_2<OfK, C>& aPoly,
                                                          const OfK& ofk,
                                                          const SsK& ssk)
{
  return arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ssk));
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_with_holes_2<OfK, C>& aPoly,
                                                          const OfK& ofk)
{
  return create_interior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, ofk, ofk);
}

template<class FT, class OfK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OfK, C> > >
inline
create_interior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_with_holes_2<OfK, C>& aPoly)
{
  return create_interior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, OfK());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

template<class FT, class OfK, class SsK, class C>
std::vector<boost::shared_ptr<Polygon_with_holes_2<OfK, C> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_2<OfK, C>& aPoly,
                                                          const OfK& ofk,
                                                          const SsK& ssk)
{
  return arrange_offset_polygons_2(
           create_exterior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk, ssk));
}

template<class FT, class OfK, class C>
std::vector<boost::shared_ptr<Polygon_with_holes_2<OfK, C> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_2<OfK, C>& aPoly,
                                                          const OfK& ofk)
{
  return create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, ofk, ofk);
}

template<class FT, class OfK, class C>
std::vector<boost::shared_ptr<Polygon_with_holes_2<OfK, C> > >
inline
create_exterior_skeleton_and_offset_polygons_with_holes_2(const FT& aOffset,
                                                          const Polygon_2<OfK, C>& aPoly)
{
  return create_exterior_skeleton_and_offset_polygons_with_holes_2(aOffset, aPoly, OfK());
}

} // end namespace CGAL

#endif
