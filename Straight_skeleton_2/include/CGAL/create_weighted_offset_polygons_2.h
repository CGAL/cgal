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
#include <CGAL/create_weighted_straight_skeleton_2.h>
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

#include <boost/range/value_type.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

namespace CGAL {

namespace CGAL_SS_i {

// ==== WARNING ====
// There is currently no way to recover simply-connectedness (see EnforceSimpleConnectedness())
// for faces that have fictitious vertices. Non-simply-connected faces can be created
// with skeletons of weighted polygons with holes.
//
// As such, you should either:
// - Use this function without holes (aHolesBegin == aHolesEnd)
// - Use this function with holes, but with weights set up such that there can be no non-simply connected
//   skeleton faces with fictitious vertices. This is for example the case when calling exterior skeleton:
//   the weights of the frame are set up to be very large and guarantee that no non-simply connected face
//   can appear.
//
// see also tags @partial_wsls_pwh

// If you are using this with some holes, you should know what you are doing
template <typename FT, typename PointIterator, typename HoleIterator,
          typename WeightIterator, typename HoleWeightsIterator,
          typename K>
std::shared_ptr< Straight_skeleton_2<K> >
create_partial_interior_weighted_straight_skeleton_2 ( const FT& aMaxTime,
                                                      PointIterator aOuterContour_VerticesBegin,
                                                      PointIterator aOuterContour_VerticesEnd,
                                                      HoleIterator aHolesBegin,
                                                      HoleIterator aHolesEnd,
                                                      WeightIterator aOuterContour_WeightsBegin,
                                                      WeightIterator aOuterContour_WeightsEnd,
                                                      HoleWeightsIterator aHoles_WeightsBegin,
                                                      HoleWeightsIterator aHoles_WeightsEnd,
                                                      const K& // aka 'SK'
                                                      )
{
  CGAL_precondition( aMaxTime > static_cast<FT>(0) ) ;

  typedef Straight_skeleton_2<K> Ss ;
  typedef Straight_skeleton_builder_traits_2<K> SsBuilderTraits;
  typedef Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

  typedef typename K::FT KFT ;

  typedef typename std::iterator_traits<PointIterator>::value_type InputPoint ;
  typedef typename Kernel_traits<InputPoint>::Kernel InputKernel ;
  typedef typename InputKernel::FT InputFT ;

  CGAL_precondition(std::distance(aOuterContour_VerticesBegin, aOuterContour_VerticesEnd) == std::distance(aOuterContour_WeightsBegin, aOuterContour_WeightsEnd));
  CGAL_precondition(std::distance(aHolesBegin, aHolesEnd) == std::distance(aHoles_WeightsBegin, aHoles_WeightsEnd));

  Cartesian_converter<InputKernel, K> conv;
  NT_converter<InputFT, KFT> wconv;

  InputFT lMaxTime = aMaxTime;
  std::optional<KFT> lOptMaxTime(conv(lMaxTime)) ;

  SsBuilder ssb( lOptMaxTime ) ;

  ssb.enter_contour ( aOuterContour_VerticesBegin, aOuterContour_VerticesEnd, conv ) ;
  ssb.enter_contour_weights ( aOuterContour_WeightsBegin, aOuterContour_WeightsEnd, wconv ) ;

  for(HoleIterator hi = aHolesBegin; hi != aHolesEnd && aHoles_WeightsBegin != aHoles_WeightsEnd; ++ hi, ++aHoles_WeightsBegin)
  {
    CGAL_precondition(std::distance(CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi)) ==
                      std::distance(aHoles_WeightsBegin->begin(), aHoles_WeightsBegin->end()));
    ssb.enter_contour(CGAL_SS_i::vertices_begin(*hi), CGAL_SS_i::vertices_end(*hi), conv);
    ssb.enter_contour_weights(aHoles_WeightsBegin->begin(), aHoles_WeightsBegin->end(), wconv);
  }

  return ssb.construct_skeleton();
}

template <typename FT, typename PointIterator, typename WeightIterator, typename K>
std::shared_ptr< Straight_skeleton_2<K> >
create_partial_exterior_weighted_straight_skeleton_2(const FT& aMaxOffset,
                                                     PointIterator aVerticesBegin,
                                                     PointIterator aVerticesEnd,
                                                     WeightIterator aWeightsBegin,
                                                     WeightIterator aWeightsEnd,
                                                     const K& k // aka 'SK'
                                                     )
{
  CGAL_precondition(aMaxOffset > 0);
  CGAL_precondition(std::distance(aWeightsBegin, aWeightsEnd) == std::distance(aVerticesBegin, aVerticesEnd));

  typedef typename std::iterator_traits<PointIterator>::value_type   Point_2;
  typedef typename Kernel_traits<Point_2>::Kernel                    IK;
  typedef typename IK::FT                                            IFT;

  std::shared_ptr<Straight_skeleton_2<K> > rSkeleton;

  // That's because we might not have FT == IK::FT (e.g. `double` and `Core`)
  // Note that we can also have IK != K (e.g. `Simple_cartesian<Core>` and `EPICK`)
  IFT lOffset = aMaxOffset;

  // @todo This likely should be done in the kernel K rather than the input kernel (i.e. the same
  // converter stuff that is done in `create_partial_exterior_straight_skeleton_2`?).
  std::optional<IFT> margin = compute_outer_frame_margin(aVerticesBegin,
                                                           aVerticesEnd,
                                                           aWeightsBegin,
                                                           aWeightsEnd,
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
    const FT frame_weight = FT(10) * *(std::max_element(aWeightsBegin, aWeightsEnd));
    CGAL_STSKEL_BUILDER_TRACE(4, "Frame weight = " << frame_weight);

    std::vector<FT> lFrameWeights(4, frame_weight);
    std::vector<std::vector<FT> > lHoleWeights;
    lHoleWeights.emplace_back(aWeightsBegin, aWeightsEnd);

    // If w[0] pointed to v_0, then when we reverse the polygon, the last polygon is pointing to v_{n-1}
    // but it is the edge v_0 v_{n-1}, which has the weight w_0.
    std::reverse(lHoleWeights[0].begin(), lHoleWeights[0].end());
    std::rotate(lHoleWeights[0].rbegin(), lHoleWeights[0].rbegin()+1, lHoleWeights[0].rend());

    // weights ensure that we cannot create a non-simply connected face with a frame halfedge
    rSkeleton = create_partial_interior_weighted_straight_skeleton_2(aMaxOffset,
                                                                     frame, frame+4,
                                                                     holes.begin(), holes.end(),
                                                                     lFrameWeights.begin(), lFrameWeights.end(),
                                                                     lHoleWeights.begin(), lHoleWeights.end(),
                                                                     k ) ;
  }

  return rSkeleton ;
}

} // namespace CGAL_SS_i

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

template<class FT, class APolygon, class HoleIterator, class Weights, class HoleWeightsIterator,
         class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< std::shared_ptr<OutPolygon> >
inline
create_interior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aOuterBoundary,
                                                        HoleIterator aHolesBegin,
                                                        HoleIterator aHolesEnd,
                                                        const Weights& aWeights,
                                                        HoleWeightsIterator aHoles_WeightsBegin,
                                                        HoleWeightsIterator aHoles_WeightsEnd,
                                                        const OfK& ofk,
                                                        const SsK& ssk)
{
  if(aHolesBegin == aHolesEnd) // see @partial_wsls_pwh
  {
    return create_offset_polygons_2<OutPolygon>(
          aOffset,
          CGAL_SS_i::dereference(
            CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
              aOffset,
              CGAL_SS_i::vertices_begin(aOuterBoundary),
              CGAL_SS_i::vertices_end  (aOuterBoundary),
              aHolesBegin,
              aHolesEnd,
              std::begin(aWeights),
              std::end(aWeights),
              aHoles_WeightsBegin,
              aHoles_WeightsEnd,
              ssk)),
          ofk);
  }
  else
  {
    return create_offset_polygons_2<OutPolygon>(
          aOffset,
          CGAL_SS_i::dereference(
            CGAL::create_interior_weighted_straight_skeleton_2(
              CGAL_SS_i::vertices_begin(aOuterBoundary),
              CGAL_SS_i::vertices_end  (aOuterBoundary),
              aHolesBegin,
              aHolesEnd,
              std::begin(aWeights),
              std::end(aWeights),
              aHoles_WeightsBegin,
              aHoles_WeightsEnd,
              ssk)),
          ofk);
  }
}

template<class FT, class APolygon, class HoleIterator, class Weights, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< std::shared_ptr<OutPolygon> >
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
std::vector< std::shared_ptr<OutPolygon> >
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
std::vector<std::shared_ptr<OutPolygon> >
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
std::vector<std::shared_ptr<OutPolygon> >
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
std::vector< std::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights,
                                                        const OfK& ofk,
                                                        const SsK& ssk,
                                                        std::enable_if_t<
                                                          ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>* = nullptr)
{
  return create_offset_polygons_2<OutPolygon>(
           aOffset,
           CGAL_SS_i::dereference(
             CGAL_SS_i::create_partial_exterior_weighted_straight_skeleton_2(
               aOffset,
               CGAL_SS_i::vertices_begin(aPoly),
               CGAL_SS_i::vertices_end  (aPoly),
               aWeights[0].begin(),
               aWeights[0].end(),
               ssk)),
           ofk);
}

// Overloads common to both polygons with and without holes, a simple polygons is returned in any case
template<class FT, class APolygon, class Weights, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< std::shared_ptr<OutPolygon> >
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
std::vector< std::shared_ptr<OutPolygon> >
inline
create_exterior_weighted_skeleton_and_offset_polygons_2(const FT& aOffset,
                                                        const APolygon& aPoly,
                                                        const Weights& aWeights)
{
  return create_exterior_weighted_skeleton_and_offset_polygons_2(aOffset, aPoly, aWeights,
                                                                 Exact_predicates_inexact_constructions_kernel());
}

} // namespace CGAL

#endif // CGAL_CREATE_WEIGHTED_OFFSET_POLYGONS_2_H
