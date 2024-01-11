// Copyright (c) 2019-2022 Google LLC (USA).
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
#ifndef CGAL_ALPHA_WRAP_3_H
#define CGAL_ALPHA_WRAP_3_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_3.h>

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
* \ingroup AW3_free_functions_grp
*
* \brief computes a watertight, 2-manifold, and intersection-free triangulated surface mesh
* that strictly contains an input triangle soup.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw3_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange a model of `Range` whose value type is the point type
* \tparam FaceRange a model of `RandomAccessContainer` whose value type is a model of `RandomAccessContainer` whose value type is an integral type
* \tparam OutputMesh model of `MutableFaceGraph`.
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \tparam OutputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param points the input points
* \param faces the input faces, with each element of the range being a range of indices corresponding to points in `points`
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output surface mesh
* \param in_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the point set `points`}
*     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
*                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
*     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
*   \cgalParamNEnd
*
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
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `alpha_wrap`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<OutputMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, alpha_wrap)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `OutputMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre The elements of `faces` are triangles.
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange, typename FaceRange, typename OutputMesh,
          typename InputNamedParameters, typename OutputNamedParameters>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap,
                  const InputNamedParameters& in_np,
                  const OutputNamedParameters& out_np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, InputNamedParameters>;
  using Geom_traits = typename NP_helper::Geom_traits;
  using Oracle = Alpha_wraps_3::internal::Triangle_soup_oracle<Geom_traits>;
  using AW3 = Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_triangle_soup(points, faces, in_np);
  AW3 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

// Convenience overloads
template <typename PointRange, typename FaceRange, typename OutputMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_3(points, faces, alpha, offset, alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename PointRange, typename FaceRange, typename OutputMesh>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap)
{
  return alpha_wrap_3(points, faces, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename PointRange, typename FaceRange, typename OutputMesh,
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputMesh& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& in_np,
                  const CGAL::Named_function_parameters<T_O, Tag_O, Base_O>& out_np,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_3(points, faces, alpha, alpha / 30., alpha_wrap, in_np, out_np);
}

template <typename PointRange, typename FaceRange, typename OutputMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputMesh& alpha_wrap,
                  const CGAL_NP_CLASS& in_np,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_3(points, faces, alpha, alpha / 30., alpha_wrap, in_np,
                      CGAL::parameters::default_values());
}

template <typename PointRange, typename FaceRange, typename OutputMesh>
void alpha_wrap_3(const PointRange& points,
                  const FaceRange& faces,
                  const double alpha,
                  OutputMesh& alpha_wrap,
                  std::enable_if_t<boost::has_range_const_iterator<FaceRange>::value>* = nullptr)
{
  return alpha_wrap_3(points, faces, alpha, alpha / 30., alpha_wrap,
                      CGAL::parameters::default_values(), CGAL::parameters::default_values());
}

// /////////////////////////////////////////////////////////////////////////////////////////////////
// WITH A TRIANGLE MESH ----------------------------------------------------------------------------
// /////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW3_free_functions_grp
*
* \brief computes a watertight, 2-manifold, and intersection-free triangulated surface mesh
* that strictly contains an input triangle mesh.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw3_parameters for a detailed breakdown of their influence.
*
* \tparam TriangleMesh model of `FaceListGraph`.
* \tparam OutputMesh model of `MutableFaceGraph`.
* \tparam InputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
* \tparam OutputNamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param tmesh a triangle mesh
* \param alpha the value of the parameter `alpha`
* \param offset the value of the parameter `offset`
* \param alpha_wrap the output surface mesh
* \param in_np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
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
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `alpha_wrap`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<OutputMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, alpha_wrap)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `OutputMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `tmesh` is a triangle mesh.
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename TriangleMesh, typename OutputMesh,
          typename InputNamedParameters, typename OutputNamedParameters>
void alpha_wrap_3(const TriangleMesh& tmesh,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap,
                  const InputNamedParameters& in_np,
                  const OutputNamedParameters& out_np
#ifndef DOXYGEN_RUNNING
                  , std::enable_if_t<! boost::has_range_const_iterator<TriangleMesh>::value>* = nullptr
#endif
                  )
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, InputNamedParameters>::type;
  using Oracle = Alpha_wraps_3::internal::Triangle_mesh_oracle<Geom_traits>;
  using AW3 = Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(alpha, gt);
  oracle.add_triangle_mesh(tmesh, in_np);
  AW3 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

// The convenience overloads are the same for triangle mesh & point set

// /////////////////////////////////////////////////////////////////////////////////////////////////
//  WITH A POINT SET -------------------------------------------------------------------------------
// /////////////////////////////////////////////////////////////////////////////////////////////////

/*!
* \ingroup AW3_free_functions_grp
*
* \brief computes a watertight, 2-manifold, and intersection-free triangulated surface mesh
* that strictly contains an input point set.
*
* The parameters `alpha` and `offset` respectively control which features will appear in the output,
* and the distance from the input. See Section \ref aw3_parameters for a detailed breakdown of their influence.
*
* \tparam PointRange model of `Range` whose value type is a point type.
* \tparam OutputMesh model of `MutableFaceGraph`.
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
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the point range}
*     \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`}
*     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
*   \cgalParamNEnd
*
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
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `alpha_wrap`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<OutputMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, alpha_wrap)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `OutputMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \pre `alpha` and `offset` are strictly positive values.
*/
template <typename PointRange, typename OutputMesh,
#ifdef DOXYGEN_RUNNING
          typename InputNamedParameters, typename OutputNamedParameters>
#else
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
#endif
void alpha_wrap_3(const PointRange& points,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap,
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
  using Oracle = Alpha_wraps_3::internal::Point_set_oracle<Geom_traits>;
  using AW3 = Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>;

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(in_np, internal_np::geom_traits));

  Oracle oracle(gt);
  oracle.add_point_set(points, in_np);
  AW3 alpha_wrap_builder(oracle);
  alpha_wrap_builder(alpha, offset, alpha_wrap, in_np, out_np);
}

// Convenience overloads, common to both mesh and point set
template <typename Input, typename OutputMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_3(const Input& input,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_3(input, alpha, offset, alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename Input, typename OutputMesh>
void alpha_wrap_3(const Input& input,
                  const double alpha,
                  const double offset,
                  OutputMesh& alpha_wrap)
{
  return alpha_wrap_3(input, alpha, offset, alpha_wrap, CGAL::parameters::default_values());
}

// without offset
template <typename Input, typename OutputMesh,
          typename T_I, typename Tag_I, typename Base_I,
          typename T_O, typename Tag_O, typename Base_O>
void alpha_wrap_3(const Input& input,
                  const double alpha,
                  OutputMesh& alpha_wrap,
                  const CGAL::Named_function_parameters<T_I, Tag_I, Base_I>& in_np,
                  const CGAL::Named_function_parameters<T_O, Tag_O, Base_O>& out_np)
{
  return alpha_wrap_3(input, alpha, alpha / 30., alpha_wrap, in_np, out_np);
}

template <typename Input, typename OutputMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
void alpha_wrap_3(const Input& input,
                  const double alpha,
                  OutputMesh& alpha_wrap,
                  const CGAL_NP_CLASS& in_np)
{
  return alpha_wrap_3(input, alpha, alpha / 30., alpha_wrap, in_np, CGAL::parameters::default_values());
}

template <typename Input, typename OutputMesh>
void alpha_wrap_3(const Input& input,
                  const double alpha,
                  OutputMesh& alpha_wrap)
{
  return alpha_wrap_3(input, alpha, alpha / 30., alpha_wrap,
                      CGAL::parameters::default_values(), CGAL::parameters::default_values());
}

} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_H
