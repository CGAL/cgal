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
#ifndef CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_2_H
#define CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/compute_outer_frame_margin.h>
#include <CGAL/Polygon_offset_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_2/Polygon_iterators.h>

#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
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

template<class FT, class PointIterator, class HoleIterator, class Weights, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_partial_interior_weighted_straight_skeleton_2 ( FT const&     aMaxTime
                                                     , PointIterator aOuterContour_VerticesBegin
                                                     , PointIterator aOuterContour_VerticesEnd
                                                     , HoleIterator  aHolesBegin
                                                     , HoleIterator  aHolesEnd
                                                     , Weights const& aWeights
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

  ssb.enter_weights( aWeights ) ;

  return ssb.construct_skeleton();
}

template<class FT, class PointIterator, class Weights, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_partial_exterior_weighted_straight_skeleton_2 ( FT const&      aMaxOffset
                                                     , PointIterator  aVerticesBegin
                                                     , PointIterator  aVerticesEnd
                                                     , Weights const& aWeights
                                                     , K const&       k // aka 'SK'
                                                    )
{
  CGAL_precondition( aMaxOffset > 0 ) ;
  CGAL_precondition(aWeights.size() == 1); // single (external) contour
  CGAL_precondition(aWeights[0].size() == std::distance(aVerticesBegin, aVerticesEnd));

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
                                                           aWeights[0],
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

    CGAL_STSKEL_BUILDER_TRACE(2, "Frame:\n" << frame[0] << "\n" << frame[1] << "\n" << frame[2] << "\n" << frame[3]);

    typedef std::vector<Point_2> Hole ;

    Hole lPoly(aVerticesBegin, aVerticesEnd);
    std::reverse(lPoly.begin(), lPoly.end());

    std::vector<Hole> holes ;
    holes.push_back(lPoly) ;

    // put a weight large enough such that frame edges are not relevant
    const FT frame_weight = *(std::max_element(aWeights[0].begin(), aWeights[0].end()));
    CGAL_STSKEL_BUILDER_TRACE(4, "Frame weight = " << frame_weight);

    Weights lWeights = { std::vector<FT>(4, frame_weight),
                         std::vector<FT>(std::rbegin(aWeights[0]), std::rend(aWeights[0])) };

    // If w[0] pointed to v_0, then when we reverse the polygon, the last polygon is pointing to v_{n-1}
    // but it is the edge v_0 v_{n-1}, which has the weight w_0.
    std::rotate(lWeights[1].rbegin(), lWeights[1].rbegin()+1, lWeights[1].rend());

    rSkeleton = create_partial_interior_weighted_straight_skeleton_2(aMaxOffset,
                                                                     frame, frame+4,
                                                                     holes.begin(), holes.end(),
                                                                     lWeights,
                                                                     k ) ;
  }

  return rSkeleton ;
}

//
// Kernel != Skeleton::kernel. The skeleton is converted to Straight_skeleton_2<Kernel>
//
template<class OutPolygon, class FT, class Skeleton, class Weights, class K>
std::vector< boost::shared_ptr<OutPolygon> >
create_weighted_offset_polygons_2 ( FT const& aOffset,
                                    Skeleton const& aSs,
                                    Weights const&, // weights are already in the HDS of the skeleton
                                    K const& ,
                                    Tag_false /*not the same kernel*/ )
{
  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ;
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;

  typedef Straight_skeleton_2<K> OfSkeleton ;

  typedef Polygon_offset_builder_traits_2<K>                                  OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<OfSkeleton,OffsetBuilderTraits,OutPolygon> OffsetBuilder;

  boost::shared_ptr<OfSkeleton> lConvertedSs = convert_straight_skeleton_2<OfSkeleton>(aSs);
  OffsetBuilder ob( *lConvertedSs );

  OutPolygonPtrVector rR ;
  ob.construct_offset_contours(aOffset, std::back_inserter(rR) ) ;

  return rR ;
}

//
// Kernel == Skeleton::kernel, no conversion
//
template<class OutPolygon, class FT, class Skeleton, class Weights, class K>
std::vector< boost::shared_ptr<OutPolygon> >
create_weighted_offset_polygons_2 ( FT const& aOffset,
                                    Skeleton const& aSs,
                                    Weights const&, // weights are already in the HDS of the skeleton
                                    K const& /*k*/,
                                    Tag_true /*same kernel*/ )
{
  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ;
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;

  typedef Polygon_offset_builder_traits_2<K>                                OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<Skeleton,OffsetBuilderTraits,OutPolygon> OffsetBuilder;

  OffsetBuilder ob(aSs);

  OutPolygonPtrVector rR ;
  ob.construct_offset_contours(aOffset, std::back_inserter(rR) ) ;

  return rR ;
}

} // namespace CGAL_SS_i

template<class OutPolygon, class FT, class Skeleton, class Weights, class K>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_weighted_offset_polygons_2(const FT& aOffset,
                                  const Skeleton& aSs,
                                  const Weights& aWeights,
                                  const K& k)
{
  typename CGAL_SS_i::Is_same_type<K, typename Skeleton::Traits>::type same_kernel;
  return CGAL_SS_i::create_weighted_offset_polygons_2<OutPolygon>(aOffset, aSs, aWeights, k, same_kernel);
}

template<class Polygon = Polygon_2<Exact_predicates_inexact_constructions_kernel>,
         class FT, class Skeleton, class Weights>
std::vector< boost::shared_ptr<Polygon> >
inline
create_weighted_offset_polygons_2(const FT& aOffset,
                                  const Skeleton& aSs,
                                  const Weights& aWeights)
{
  return create_weighted_offset_polygons_2<Polygon>(aOffset, aSs, aWeights, Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

template<class FT, class APolygon, class HoleIterator, class Weights, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aOuterBoundary,
                                                        HoleIterator aHolesBegin,
                                                        HoleIterator aHolesEnd,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk)
{
  return create_weighted_offset_polygons_2<OutPolygon>(
           aOffset,
           CGAL_SS_i::dereference(
             CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
               aOffset,
               CGAL_SS_i::vertices_begin(aOuterBoundary),
               CGAL_SS_i::vertices_end  (aOuterBoundary),
               aWeights,
               aHolesBegin,
               aHolesEnd,
               ssk)),
           aWeights,
           ofk);
}

template<class FT, class APolygon, class HoleIterator, class Weights, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aOuterBoundary,
                                               HoleIterator aHolesBegin,
                                               HoleIterator aHolesEnd,
                                                const Weights& aWeights,
                                               const OfK& ofk)
{
  return create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aOuterBoundary,
                                                                 aHolesBegin, aHolesEnd,
                                                                 aWeights,
                                                                 ofk,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

// Overload where Polygon actually is a simple polygon (no holes)
template<class FT, class APolygon, class Weights, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk,
                                                        std::enable_if_t<
                                                          ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>* = nullptr)
{
  std::vector<APolygon> no_holes;
  return create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly,
                                                                 no_holes.begin(), no_holes.end(),
                                                                 aWeights,
                                                                 ofk, ssk);
}

// Overloads common to both polygons with and without holes, a simple polygon is returned in any case
template<class FT, class APolygon, class Weights, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk)
{
  return create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights, ofk,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class APolygon, class Weights,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<
                              APolygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights)
{
  return create_interior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_exterior_skeleton_and_offset_polygons_2 (no sorting of the result) */

// Overload where Polygon actually is a simple polygon (no holes)
template<class FT, class APolygon, class Weights, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk,
                                                        std::enable_if_t<
                                                          ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>* = nullptr)
{
  return create_weighted_offset_polygons_2<OutPolygon>(
           aOffset,
           CGAL_SS_i::dereference(
             CGAL_SS_i::create_partial_exterior_weighted_straight_skeleton_2(
               aOffset,
               CGAL_SS_i::vertices_begin(aPoly),
               CGAL_SS_i::vertices_end  (aPoly),
               aWeights,
               ssk)),
           aWeights,
           ofk);
}

// Overloads common to both polygons with and without holes, a simple polygons is returned in any case
template<class FT, class APolygon, class Weights, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk)
{
  return create_exterior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights, ofk,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class APolygon, class Weights,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<
                              APolygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights)
{
  return create_exterior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

} // end namespace CGAL

#endif
