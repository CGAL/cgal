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

#include <CGAL/create_straight_skeleton_2.h>
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

template<class U, class V>
struct Is_same_type { typedef Tag_false type ; } ;

template<class U>
struct Is_same_type<U,U> { typedef Tag_true type ; } ;

// The return type of create_interior/exterior_skeleton_and_offset_polygons_2:
// - if polygon input is a model of 'GeneralPolygonWithHoles_2', the return type
//   should be the internal (hole-less) polygon type GeneralPolygonWithHoles_2::General_polygon_2
// - if polygon input is just a sequence container of points (e.g. Polygon_2), then the same type
//   is expected in output
template <typename Polygon, typename OfK,
          bool has_holes = CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>
struct Default_return_polygon_type // Polygon type supports holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<
                                        typename Polygon::General_polygon_2>::type>::Kernel,
                                      OfK>::value,
                                    typename Polygon::General_polygon_2, // correct kernel
                                    CGAL::Polygon_2<OfK> /*incorrect kernel*/ >::type type;
};

template <typename Polygon, typename OfK>
struct Default_return_polygon_type<Polygon, OfK, false> // Polygon type does NOT support holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<Polygon>::type>::Kernel,
                                      OfK>::value,
                                    Polygon, // correct kernel
                                    CGAL::Polygon_2<OfK> /*incorrect kernel*/ >::type type;
};

// The return type of create_interior/exterior_skeleton_and_offset_polygons_with_holes_2:
// - if polygon input is a model of 'GeneralPolygonWithHoles_2', the return type should be the same
// - if polygon input is just a sequence container of points (e.g. Polygon_2), then use
//   General_polygon_with_holes_2<Polygon>
template <typename Polygon, typename OfK,
          bool has_holes = CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>
struct Default_return_polygon_with_holes_type // Polygon type supports holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<
                                        typename Polygon::General_polygon_2>::type>::Kernel,
                                      OfK>::value,
                                    Polygon, // correct kernel
                                    CGAL::Polygon_with_holes_2<OfK> /*incorrect kernel*/ >::type type;
};

template <typename Polygon, typename OfK>
struct Default_return_polygon_with_holes_type<Polygon, OfK, false> // Polygon type does NOT support holes
{
  // Maybe on paper the `conditional<true>` should be `General_polygon_with_holes_2<Polygon>`...
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<Polygon>::type>::Kernel,
                                      OfK>::value,
                                    CGAL::Polygon_with_holes_2<OfK>, // correct kernel but no holes
                                    CGAL::Polygon_with_holes_2<OfK> /*incorrect kernel*/ >::type type;
};

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
    double lm = CGAL::to_double(*margin);
    Bbox_2 bbox = bbox_2(aVerticesBegin, aVerticesEnd);

    FT fxmin = bbox.xmin() - lm ;
    FT fxmax = bbox.xmax() + lm ;
    FT fymin = bbox.ymin() - lm ;
    FT fymax = bbox.ymax() + lm ;

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

    rSkeleton = create_partial_interior_straight_skeleton_2(aMaxOffset,frame, frame+4, holes.begin(), holes.end(), k ) ;
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
  typedef boost::shared_ptr<OutPolygon> OutPolygonPtr ;
  typedef std::vector<OutPolygonPtr>    OutPolygonPtrVector ;

  typedef Polygon_offset_builder_traits_2<K>                                OffsetBuilderTraits;
  typedef Polygon_offset_builder_2<Skeleton,OffsetBuilderTraits,OutPolygon> OffsetBuilder;

  OutPolygonPtrVector rR ;

  OffsetBuilder ob(aSs);
  ob.construct_offset_contours(aOffset, std::back_inserter(rR) ) ;

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

template<class OutPolygon, class FT, class Skeleton, class K>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_offset_polygons_2(const FT& aOffset,
                         const Skeleton& aSs,
                         const K& k)
{
  typename CGAL_SS_i::Is_same_type<K, typename Skeleton::Traits>::type same_kernel;
  return CGAL_SS_i::create_offset_polygons_2<OutPolygon>(aOffset, aSs, k, same_kernel);
}

template<class Polygon = Polygon_2<Exact_predicates_inexact_constructions_kernel>,
         class FT, class Skeleton>
std::vector< boost::shared_ptr<Polygon> >
inline
create_offset_polygons_2(const FT& aOffset,
                         const Skeleton& aSs)
{
  return create_offset_polygons_2<Polygon>(aOffset, aSs, Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// INTERIOR

template<class FT, class APolygon, class HoleIterator, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aOuterBoundary,
                                               HoleIterator aHolesBegin,
                                               HoleIterator aHolesEnd,
                                               const OfK& ofk,
                                               const SsK& ssk)
{
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

template<class FT, class APolygon, class HoleIterator, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aOuterBoundary,
                                               HoleIterator aHolesBegin,
                                               HoleIterator aHolesEnd,
                                               const OfK& ofk)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aOuterBoundary,
                                                        aHolesBegin, aHolesEnd,
                                                        ofk,
                                                        Exact_predicates_inexact_constructions_kernel());
}

// Overload where Polygon actually is a simple polygon (no holes)
template<class FT, class APolygon, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk,
                                               const SsK& ssk,
                                               typename std::enable_if<
                                                 ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>::type* = nullptr)
{
  std::vector<APolygon> no_holes;
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly,
                                                        no_holes.begin(), no_holes.end(),
                                                        ofk, ssk);
}

// Overloads common to both polygons with and without holes, a simple polygon is returned in any case
template<class FT, class APolygon, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk,
                                                        Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class APolygon,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<
                              APolygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector<boost::shared_ptr<OutPolygon> >
inline
create_interior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly)
{
  return create_interior_skeleton_and_offset_polygons_2(aOffset, aPoly,
                                                        Exact_predicates_inexact_constructions_kernel());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// EXTERIOR

/*! create_exterior_skeleton_and_offset_polygons_2 (no sorting of the result) */

// Overload where Polygon actually is a simple polygon (no holes)
template<class FT, class APolygon, class OfK, class SsK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk,
                                               const SsK& ssk,
                                               typename std::enable_if<
                                                 ! CGAL_SS_i::has_Hole_const_iterator<APolygon>::value>::type* = nullptr)
{
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

// Overloads common to both polygons with and without holes, a simple polygon is returned in any case
template<class FT, class APolygon, class OfK,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<APolygon, OfK>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly,
                                               const OfK& ofk)
{
  return create_exterior_skeleton_and_offset_polygons_2(aOffset, aPoly, ofk,
                                                        Exact_predicates_inexact_constructions_kernel());
}

template<class FT, class APolygon,
         class OutPolygon = typename CGAL_SS_i::Default_return_polygon_type<
                              APolygon, Exact_predicates_inexact_constructions_kernel>::type>
std::vector< boost::shared_ptr<OutPolygon> >
inline
create_exterior_skeleton_and_offset_polygons_2(const FT& aOffset,
                                               const APolygon& aPoly)
{
  return create_exterior_skeleton_and_offset_polygons_2(aOffset, aPoly,
                                                        Exact_predicates_inexact_constructions_kernel());
}

} // end namespace CGAL

#endif
