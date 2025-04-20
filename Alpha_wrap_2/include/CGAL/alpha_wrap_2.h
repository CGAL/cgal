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

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_2.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <boost/range/has_range_iterator.hpp>

#include <type_traits>

namespace CGAL {

// /////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A TRIANGLE SOUP ---------------------------------------------------------------------------
// /////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a range of simple polygons that strictly contain an input triangle soup.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange a model of `Range` whose value type is the point type
* \tparam FaceRange a model of `RandomAccessContainer` whose value type is a model of `RandomAccessContainer` whose value type is an integral type
* \tparam OutputPolygons model of `BackInsertionSequence` whose value is a model of `GeneralPolygonWithHoles_2`.
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \tparam OutputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param points the input points
* \param faces the input faces, with each element of the range being a range of indices corresponding to points in `points`
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output polygons
* \param in_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
* \param out_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \pre The elements of `faces` are triangles.
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange, typename FaceRange, typename OutputPolygons,
          typename InputNamedParameters, typename OutputNamedParameters>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap,
                  const InputNamedParameters& in_np,
                  const OutputNamedParameters& out_np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, InputNamedParameters>;
  using Geom_traits = typename NP_helper::Geom_traits;
  using Oracle = Alpha_wraps_2::internal::Triangle_soup_oracle<Geom_traits>;
  using AW2 = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_triangle_soup(points, faces, in_np);
  AW2 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

// Convenience overloads
template <typename PointRange, typename FaceRange, typename OutputPolygons,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_2(points, faces, alpha, offset, alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename PointRange, typename FaceRange, typename OutputPolygons>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap)
{
  return alpha_wrap_2(points, faces, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename PointRange, typename FaceRange, typename OutputPolygons,
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputPolygons& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& in_np,
                  const CGAL::Named_function_parameters<T_O, Tag_O, Base_O>& out_np,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_2(points, faces, alpha, alpha / 30., alpha_wrap, in_np, out_np);
}

template <typename PointRange, typename FaceRange, typename OutputPolygons,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputPolygons& alpha_wrap,
                  const CGAL_NP_CLASS& in_np,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_2(points, faces, alpha, alpha / 30., alpha_wrap, in_np,
                      CGAL::parameters::default_values());
}

template <typename PointRange, typename FaceRange, typename OutputPolygons>
void alpha_wrap_2(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputPolygons& alpha_wrap,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_2(points, faces, alpha, alpha / 30., alpha_wrap,
                      CGAL::parameters::default_values(), CGAL::parameters::default_values());
}

#if 0

// /////////////////////////////////////////////////////////////////////////////////////////////////
// WITH POLYGONS -----------------------------------------------------------------------------------
// /////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a range of simple polygons that strictly contain a range of polygons.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam InputPolygons model of `Range` whose value type is a model of `GeneralPolygonWithHoles_2`.
* \tparam OutputPolygons model of `BackInsertionSequence` whose value type is a model of `GeneralPolygonWithHoles_2`.
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \tparam OutputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param polygons a range of polygons
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output surface mesh
* \param in_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
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
template <typename InputPolygons, typename OutputPolygons,
          typename InputNamedParameters, typename OutputNamedParameters>
void alpha_wrap_2(const InputPolygons& polygons,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap,
                  const InputNamedParameters& in_np,
                  const OutputNamedParameters& out_np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using Geom_traits = typename boost::range_value<InputPolygons>::type::Geom_traits;
  using Oracle = Alpha_wraps_2::internal::Polygons_oracle<Geom_traits>;
  using AW2 = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_polygons(polygons, in_np);
  AW2 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

#endif

// The convenience overloads are the same for polygons & point set

// /////////////////////////////////////////////////////////////////////////////////////////////////
//  WITH A POINT SET -------------------------------------------------------------------------------
// /////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW2_free_functions_grp
*
* \brief computes a watertight, 2-manifold, and intersection-free triangulated surface mesh
* that strictly contains an input point set.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw2_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange model of `Range` whose value type is a point type.
* \tparam OutputPolygons model of `BackInsertionSequence` whose value is a model of `GeneralPolygonWithHoles_2`.
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \tparam OutputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param points the input points
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output surface mesh
* \param in_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{<ul><li>The geometric traits class must be compatible with the point type.</li>
*                         <li>The geometric traits should use a floating point number type (see \ref aw3_interface).</li></ul>}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \param out_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange, typename OutputPolygons,
#ifdef DOXYGEN_RUNNING
          typename InputNamedParameters, typename OutputNamedParameters>
#else
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
#endif
void alpha_wrap_2(const PointRange& points,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap,
#ifdef DOXYGEN_RUNNING
                  const InputNamedParameters& in_np,
                  const OutputNamedParameters& out_np
#else
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& in_np,
                  const CGAL::Named_function_parameters<T_O, Tag_O, Base_O>& out_np,
                  std::enable_if_t<boost::has_range_const_iterator<PointRange>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using InputNamedParameters = CGAL::Named_function_parameters<T_I, Tag_I, Base_I>;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, InputNamedParameters>;
  using Geom_traits = typename NP_helper::Geom_traits;
  using Oracle = Alpha_wraps_2::internal::Point_set_oracle<Geom_traits>;
  using AW2 = Alpha_wraps_2::internal::Alpha_wrapper_2<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(gt);
  oracle.add_point_set(points, in_np);
  AW2 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

// Convenience overloads, common to both polygons and point set
template <typename Input, typename OutputPolygons,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_2(input, alpha, offset, alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename Input, typename OutputPolygons>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  const double offset,
                  OutputPolygons& alpha_wrap)
{
  return alpha_wrap_2(input, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename Input, typename OutputPolygons,
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  OutputPolygons& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& in_np,
                  const CGAL::Named_function_parameters<T_O, Tag_O, Base_O>& out_np)
{
  return alpha_wrap_2(input, alpha, alpha / 30., alpha_wrap, in_np, out_np);
}

template <typename Input, typename OutputPolygons, typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  OutputPolygons& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_2(input, alpha, alpha / 30., alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename Input, typename OutputPolygons>
void alpha_wrap_2(const Input& input,
                  const double alpha,
                  OutputPolygons& alpha_wrap)
{
  return alpha_wrap_2(input, alpha, alpha / 30., alpha_wrap,
                      CGAL::parameters::default_values(), CGAL::parameters::default_values());
}

} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_H
