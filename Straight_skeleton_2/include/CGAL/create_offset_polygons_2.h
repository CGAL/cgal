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
#ifndef CGAL_CREATE_OFFSET_POLYGONS_2_H
#define CGAL_CREATE_OFFSET_POLYGONS_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_2/Polygon_iterators.h>

#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Default.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/tags.h>

#include <boost/optional/optional.hpp>
#include <boost/range/value_type.hpp>
#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

namespace CGAL {

namespace CGAL_SS_i {

template<class FT, class PointIterator, class HoleIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_partial_interior_straight_skeleton_2 ( FT const&     aMaxTime
                                            , PointIterator aOuterContour_VerticesBegin
                                            , PointIterator aOuterContour_VerticesEnd
                                            , HoleIterator  aHolesBegin
                                            , HoleIterator  aHolesEnd
                                            , K const& // aka 'SK'
                                            )
{
  CGAL_precondition( aMaxTime > static_cast<FT>(0) ) ;

  typedef Straight_skeleton_2<K> Ss ;
  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

  typedef typename K::FT KFT ;

  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;

  Cartesian_converter<InputKernel, K> conv ;

  typename InputKernel::FT lMaxTime = aMaxTime;
  boost::optional<KFT> lOptMaxTime(conv(lMaxTime)) ;

  SsBuilder ssb( lOptMaxTime ) ;

  ssb.enter_contour( aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, conv ) ;

  for ( HoleIterator hi = aHolesBegin ; hi != aHolesEnd ; ++ hi )
    ssb.enter_contour( CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi), conv ) ;

  return ssb.construct_skeleton();
}

template<class FT, class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_partial_exterior_straight_skeleton_2 ( FT const&      aMaxOffset
                                            , PointIterator  aVerticesBegin
                                            , PointIterator  aVerticesEnd
                                            , K const&       k // aka 'SK'
                                            )
{
  CGAL_precondition( aMaxOffset > 0 ) ;

  typedef typename std::iterator_traits<PointIterator>::value_type   Point_2;
  typedef typename Kernel_traits<Point_2>::Kernel                    IK;
  typedef typename IK::FT                                            IFT;

  boost::shared_ptr<Straight_skeleton_2<K> > rSkeleton;

  // That's because we might not have FT == IK::FT (e.g. `double` and `Core`)
  // Note that we can also have IK != K (e.g. `Simple_cartesian<Core>` and `EPICK`)
  IFT lOffset = aMaxOffset;

  // @todo This likely should be done in the kernel K rather than the input kernel (i.e. the same
  // converter stuff that is done in `create_partial_exterior_straight_skeleton_2`?).
  boost::optional<IFT> margin = compute_outer_frame_margin(aVerticesBegin,
                                                           aVerticesEnd,
                                                           lOffset);

  if ( margin )
  {
    const IFT lm = *margin;
    const Bbox_2 bbox = bbox_2(aVerticesBegin, aVerticesEnd);

    const IFT fxmin = IFT(bbox.xmin()) - lm ;
    const IFT fxmax = IFT(bbox.xmax()) + lm ;
    const IFT fymin = IFT(bbox.ymin()) - lm ;
    const IFT fymax = IFT(bbox.ymax()) + lm ;

    Point_2 frame[4] ;
    frame[0] = Point_2(fxmin,fymin) ;
    frame[1] = Point_2(fxmax,fymin) ;
    frame[2] = Point_2(fxmax,fymax) ;
    frame[3] = Point_2(fxmin,fymax) ;

    typedef std::vector<Point_2> Hole ;

    Hole lPoly(aVerticesBegin, aVerticesEnd);
    std::reverse(lPoly.begin(), lPoly.end());

    std::vector<Hole> holes ;
    holes.push_back(lPoly) ;

    rSkeleton = create_partial_interior_straight_skeleton_2(aMaxOffset, frame, frame+4, holes.begin(), holes.end(), k ) ;
  }

  return rSkeleton ;
}

//
// Kernel != Skeleton::kernel. The skeleton is converted to Straight_skeleton_2<Kernel>
//
template<class OutPolygon, class FT, class Skeleton, class K>
std::vector< boost::shared_ptr<OutPolygon> >
create_offset_polygons_2 ( FT const& aOffset, Skeleton const& aSs, K const& , Tag_false )
{
  static_assert(!std::is_same_v<OutPolygon, CGAL::Default>);

  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ;
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;

  typedef Straight_skeleton_2<K> OfSkeleton ;

  typedef Polygon_offset_builder_traits_2<K>                                  OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<OfSkeleton,OffsetBuilderTraits,OutPolygon> OffsetBuilder;

  OutPolygonPtrVector rR ;

  boost::shared_ptr<OfSkeleton> lConvertedSs = convert_straight_skeleton_2<OfSkeleton>(aSs);
  OffsetBuilder ob( *lConvertedSs );
  ob.construct_offset_contours(aOffset, std::back_inserter(rR) ) ;

  return rR ;
}

//
// Kernel == Skeleton::kernel, no conversion
//
template<class OutPolygon, class FT, class Skeleton, class K>
std::vector< boost::shared_ptr<OutPolygon> >
create_offset_polygons_2 ( FT const& aOffset, Skeleton const& aSs, K const& /*k*/, Tag_true )
{
  static_assert(!std::is_same_v<OutPolygon, CGAL::Default>);

  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ;
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;

  typedef Polygon_offset_builder_traits_2<K>                                OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<Skeleton,OffsetBuilderTraits,OutPolygon> OffsetBuilder;

  OffsetBuilder ob(aSs);
  typename K::FT lOffset = aOffset;
  OutPolygonPtrVector rR ;
  ob.construct_offset_contours(lOffset, std::back_inserter(rR) ) ;

  return rR ;
}

// Allow failure due to invalid straight skeletons to go through the users
template<class Skeleton>
Skeleton const& dereference ( boost::shared_ptr<Skeleton> const& ss )
{
  CGAL_precondition(ss.get() != 0);
  return *ss;
}

} // namespace CGAL_SS_i

template<class OutPolygon, class FT, class Skeleton,
         class K = Exact_predicates_inexact_constructions_kernel>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_offset_polygons_2(const FT& aOffset,
                         const Skeleton& aSs,
                         const K& k = K())
{
  typename CGAL_SS_i::Is_same_type<K, typename Skeleton::Traits>::type same_kernel;
  return CGAL_SS_i::create_offset_polygons_2<OutPolygon>(aOffset, aSs, k, same_kernel);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

template <class OutPolygon_ = CGAL::Default,
          class FT, class APolygon, class HoleIterator,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector< boost::shared_ptr<CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aOuterBoundary,
                                               HoleIterator aHolesBegin,
                                               HoleIterator aHolesEnd,
                                               const OfK& ofk = OfK(),
                                               const SsK& ssk = SsK(),
                                               std::enable_if_t<CGAL::is_iterator<HoleIterator>::value>* = 0)
{
  using OutPolygon = CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK>;

  return create_offset_polygons_2<OutPolygon>(
           aOffset,
           CGAL_SS_i::dereference(
             CGAL_SS_i::create_partial_interior_straight_skeleton_2(
               aOffset,
               CGAL_SS_i::vertices_begin(aOuterBoundary),
               CGAL_SS_i::vertices_end  (aOuterBoundary),
               aHolesBegin,
               aHolesEnd,
               ssk)),
           ofk);
}

// Overload where APolygon is a simple polygon (no holes)
template <class OutPolygon_ = CGAL::Default,
          class FT, class APolygon,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector< boost::shared_ptr<CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK> > >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk = OfK(),
                                               const SsK& ssk = SsK(),
                                               std::enable_if_t<
                                                 ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>* = nullptr)
{
  using OutPolygon = CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK>;

  std::vector<APolygon> no_holes;
  return create_interior_skeleton_and_offset_polygons_2<OutPolygon>(aOffset, aPoly,
                                                                    no_holes.begin(), no_holes.end(),
                                                                    ofk, ssk);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_exterior_skeleton_and_offset_polygons_2 (no sorting of the result) */

// Overload where Polygon actually is a simple polygon (no holes)
template <class OutPolygon_ = CGAL::Default,
          class FT, class APolygon,
          class OfK = Exact_predicates_inexact_constructions_kernel,
          class SsK = Exact_predicates_inexact_constructions_kernel>
std::vector< boost::shared_ptr<CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK> > >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk = OfK(),
                                               const SsK& ssk = SsK(),
                                               std::enable_if_t<
                                                 ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>* = nullptr)
{
  using OutPolygon = CGAL_SS_i::Polygon_return_type<OutPolygon_, APolygon, OfK>;

  return create_offset_polygons_2<OutPolygon>(
           aOffset,
           CGAL_SS_i::dereference(
             CGAL_SS_i::create_partial_exterior_straight_skeleton_2(
               aOffset,
               CGAL_SS_i::vertices_begin(aPoly),
               CGAL_SS_i::vertices_end  (aPoly),
               ssk)),
           ofk);
}

} // namespace CGAL

#endif // CGAL_CREATE_OFFSET_POLYGONS_2_H
