// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_H
#define CGAL_ALPHA_WRAP_2_H

// #define CGAL_AW2_USE_SORTED_PRIORITY_QUEUE

#ifdef CGAL_AW2_DEBUG_PP
 #ifndef CGAL_AW2_DEBUG
  #define CGAL_AW2_DEBUG
  #define CGAL_AW2_DEBUG_INITIALIZATION
  #define CGAL_AW2_DEBUG_STEINER_COMPUTATION
  #define CGAL_AW2_DEBUG_QUEUE_PP
  #define CGAL_AW2_DEBUG_QUEUE
  #define CGAL_AW2_DEBUG_EDGE_STATUS
  #define CGAL_AW2_DEBUG_MANIFOLDNESS
  #define CGAL_AW2_DEBUG_TRAVERSABILITY
  #define CGAL_AW2_DEBUG_SPHERE_MARCHING

  #define CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
  #define CGAL_AW2_DEBUG_DUMP_EVERY_STEP
 #endif
#endif

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_2.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <boost/range/has_range_iterator.hpp>
#include <boost/range/value_type.hpp>

#include <type_traits>

namespace CGAL {

namespace Alpha_wraps_2 {
namespace internal {

// We must first check that Kernel_traits does not return the dummy kernel,
// otherwise writing `Kernel_traits<value_type>::type::Point_2` is a real compilation error.

#define CGAL_IS_RANGE_OF_KERNEL_OBJECT(TYPE)                                                       \
template <typename Input>                                                                          \
struct is_##TYPE##_range                                                                           \
{                                                                                                  \
  static constexpr bool value = []() -> bool {                                                     \
    if constexpr (!boost::has_range_const_iterator<Input>::value) {                                \
      return false;                                                                                \
    } else {                                                                                       \
      using raw_value_type = typename boost::range_value<Input>::type;                             \
      using value_type = CGAL::cpp20::remove_cvref_t<raw_value_type>;                              \
      using Kernel = typename CGAL::Kernel_traits<value_type>::type;                               \
      using Dummy_kernel = CGAL::internal_kernel_traits::Dummy_kernel<value_type>;                 \
      if constexpr(std::is_same_v<Kernel, Dummy_kernel>)                                           \
        return false;                                                                              \
      else                                                                                         \
        return std::is_same_v<value_type, typename Kernel::TYPE>;                                  \
    }                                                                                              \
  }();                                                                                             \
};

CGAL_IS_RANGE_OF_KERNEL_OBJECT(Point_2)
CGAL_IS_RANGE_OF_KERNEL_OBJECT(Segment_2)
CGAL_IS_RANGE_OF_KERNEL_OBJECT(Triangle_2)

#undef CGAL_IS_RANGE_OF_KERNEL_OBJECT

// multipolygon, check for a typedef 'Polygon_with_holes_2' required by `MultipolygonWithHoles_2`
template <typename Input>
struct is_MultipolygonWithHoles
{
  BOOST_MPL_HAS_XXX_TRAIT_DEF(Polygon_with_holes_2)
  static constexpr bool value = has_Polygon_with_holes_2<Input>::value;
};

// A multi-linestring is a range of ranges of points
template <typename Input>
struct is_MultiLineString
{
  static constexpr bool value = [] -> bool {
    if constexpr (boost::has_range_const_iterator<Input>::value) {
      using Outer_range = typename boost::range_value<Input>::type;
      return is_Point_2_range<Outer_range>::value;
    } else {
      return false;
    }
  }();
};

} // namespace internal
} // namespace Alpha_wraps_2

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH AN INDEXED FACE SET ------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* an input indexed face set.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange a model of `RandomAccessContainer` whose value type is a point type model of `Kernel::Point_2`
* \tparam FaceRange a model of `Range` whose value type is a model of `RandomAccessContainer`
*                   whose value type is an integral type
* \tparam MultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param points the input points
* \param faces the input faces, with each element of the range being a range of indices corresponding to points in `points`
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the point type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange, typename FaceRange,
          typename MultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, InputNamedParameters>;
  using Geom_traits = typename NP_helper::Geom_traits;
  using Oracle = Alpha_wraps_2::internal::Segment_soup_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_polygon_soup(points, faces, np);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, np);
}

// Convenience overloads

template <typename PointRange, typename FaceRange, typename MultipolygonWithHoles>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap)
{
  return alpha_wrap_2(points, faces, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename PointRange, typename FaceRange, typename MultipolygonWithHoles,
          typename T_I, typename Tag_I, typename Base_I>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  MultipolygonWithHoles& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& np,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_2(points, faces, alpha, alpha / 30., alpha_wrap, np);
}

template <typename PointRange, typename FaceRange, typename MultipolygonWithHoles>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  MultipolygonWithHoles& alpha_wrap,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_2(points, faces, alpha, alpha / 30., alpha_wrap, CGAL::parameters::default_values());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A TRIANGLE SOUP ----------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* an input triangle soup.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam TriangleRange a model of `Range` whose value type is a model of `RandomAccessContainer`
*                       of size 3 whose value type is a point type model of `Kernel::Point_2`
* \tparam MultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param triangles the input triangles
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the triangle type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename TriangleRange,
          typename MultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const TriangleRange& triangles,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<Alpha_wraps_2::internal::is_Triangle_2_range<TriangleRange>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using In_triangle_2 = typename boost::range_value<TriangleRange>::type;
  using In_K = typename CGAL::Kernel_traits<In_triangle_2>::type;

  using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                   InputNamedParameters,
                                                                   In_K>::type;

  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Out_point_2 = typename boost::range_value<Polygon_2>::type;
  using Out_K = typename CGAL::Kernel_traits<Out_point_2>::type;

  // could imagine this being just a conversion, but ask for equality for now
  static_assert(std::is_same_v<In_K, Out_K>);

  using Oracle = Alpha_wraps_2::internal::Segment_soup_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_triangles(triangles);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap);
}

// The convenience overloads are common to all ranges

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A MULTI-POLYGON ----------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* an input range of polygons.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \note This function merely inserts and wraps all edges from the outer boundaries and from the
* holes of the multipolygon. It does *not* repair the multipolygon first. Consequently, for example,
* if a polygon were to be invalid with an edge that lays outside of the region enclosed by its outer
* boundary, then the wrap would still be constructed as to contain this edge.
* In some use cases, it may be desirable to first repair invalid polygons. We refer
* to the \cgal component \link Chapter_2D_Polygon_repair Polygon Repair\endlink for various
* polygon repair functions using different strategies.
*
* \tparam InputMultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam OutputMultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param multipolygon a multipolygon
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits` and
*                       the point type of the multipolygon}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the point type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename InputMultipolygonWithHoles,
          typename OutputMultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const InputMultipolygonWithHoles& multipolygon,
                  const double alpha,
                  const double offset,
                  OutputMultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<Alpha_wraps_2::internal::is_MultipolygonWithHoles<InputMultipolygonWithHoles>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using In_polygon_with_holes_2 = typename InputMultipolygonWithHoles::Polygon_with_holes_2;
  using In_polygon_2 = typename In_polygon_with_holes_2::Polygon_2;
  using In_point_2 = typename boost::range_value<In_polygon_2>::type;
  using In_K = typename CGAL::Kernel_traits<In_point_2>::type;

  using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                   InputNamedParameters,
                                                                   In_K>::type;

  using Out_polygon_with_holes_2 = typename OutputMultipolygonWithHoles::Polygon_with_holes_2;
  using Out_polygon_2 = typename Out_polygon_with_holes_2::Polygon_2;
  using Out_point_2 = typename boost::range_value<Out_polygon_2>::type;
  using Out_K = typename CGAL::Kernel_traits<Out_point_2>::type;

  // could imagine this being just a conversion, but ask for equality for now
  static_assert(std::is_same_v<In_K, Out_K>);

  using Oracle = Alpha_wraps_2::internal::Segment_soup_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_multipolygon(multipolygon);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap);
}

// The convenience overloads are common to all ranges

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A SEGMENT SOUP -----------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* an input segment soup.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam SegmentRange a model of `Range` whose value type is a segment type model of `Kernel::Segment_2`
* \tparam MultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segments
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the segment type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename SegmentRange,
          typename MultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const SegmentRange& segments,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<Alpha_wraps_2::internal::is_Segment_2_range<SegmentRange>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using In_segment_2 = typename boost::range_value<SegmentRange>::type;
  using In_K = typename CGAL::Kernel_traits<In_segment_2>::type;

  using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                   InputNamedParameters,
                                                                   In_K>::type;

  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Out_point_2 = typename boost::range_value<Polygon_2>::type;
  using Out_K = typename CGAL::Kernel_traits<Out_point_2>::type;

  using Oracle = Alpha_wraps_2::internal::Segment_soup_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_segments(segments);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A MULTI-POLYLINE ---------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* input polylines.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam MultiLineString a model of `RandomAccessContainer` whose value type is a model of `RandomAccessContainer`
* whose value type is a point type model of `Kernel::Point_2`
* \tparam MultipolygonWithHoles a model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param multilinestring the input polylines
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the segment type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename MultiLineString,
          typename MultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const MultiLineString& multilinestring,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<Alpha_wraps_2::internal::is_MultiLineString<MultiLineString>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using Linestring = typename boost::range_value<MultiLineString>::type;
  using In_Point_2 = typename boost::range_value<Linestring>::type;
  using In_K = typename CGAL::Kernel_traits<In_Point_2>::type;

  using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                   InputNamedParameters,
                                                                   In_K>::type;

  using Polygon_with_holes_2 = typename MultipolygonWithHoles::Polygon_with_holes_2;
  using Polygon_2 = typename Polygon_with_holes_2::Polygon_2;
  using Out_point_2 = typename boost::range_value<Polygon_2>::type;
  using Out_K = typename CGAL::Kernel_traits<Out_point_2>::type;

  using Oracle = Alpha_wraps_2::internal::Segment_soup_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_multilinestring(multilinestring);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap);
}

// The convenience overloads are common to all ranges

////////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A POINT SET / MULTI POINT ------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 1-manifold, simple multipolygon that strictly contains
* an input point set.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange model of `Range` whose value type is a point type model of `Kernel::Point_2`
* \tparam MultipolygonWithHoles model of `MultipolygonWithHoles_2`
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param points the input points
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output multipolygon
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the point type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw2_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange,
          typename MultipolygonWithHoles,
          typename InputNamedParameters>
void alpha_wrap_2(const PointRange& points,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap,
                  const InputNamedParameters& np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<Alpha_wraps_2::internal::is_Point_2_range<PointRange>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, InputNamedParameters>;
  using Geom_traits = typename NP_helper::Geom_traits;
  using Oracle = Alpha_wraps_2::internal::Point_set_oracle<Geom_traits>;
  using Wrapper = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  Oracle oracle(gt);
  oracle.add_points(points, np);

  Wrapper alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, np);
}

// Convenience overloads, common to all ranges

template <typename Input, typename MultipolygonWithHoles>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  const double offset,
                  MultipolygonWithHoles& alpha_wrap)
{
  return alpha_wrap_2(input, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename Input, typename MultipolygonWithHoles,
          typename T_I, typename Tag_I, typename Base_I>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  MultipolygonWithHoles& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& np)
{
  return alpha_wrap_2(input, alpha, alpha / 30., alpha_wrap, np);
}

template <typename Input, typename MultipolygonWithHoles>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  MultipolygonWithHoles& alpha_wrap)
{
  return alpha_wrap_2(input, alpha, alpha / 30., alpha_wrap, CGAL::parameters::default_values());
}

} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_H
