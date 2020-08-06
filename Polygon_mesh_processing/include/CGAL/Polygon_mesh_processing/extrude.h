// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Maxime Gimeno


#ifndef CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H
#define CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace extrude_impl{

template<typename BorderHalfedgesRange, class PolygonMesh>
void create_strip(const BorderHalfedgesRange& input_halfedges,
                 const BorderHalfedgesRange& output_halfedges,
                  PolygonMesh& mesh)
{
  CGAL_assertion(input_halfedges.size() == output_halfedges.size());
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  for(std::size_t i = 0; i < input_halfedges.size(); ++i)
  {
    halfedge_descriptor h1 = input_halfedges[i], h2=output_halfedges[i],
        nh1 = next(h1, mesh), ph2 = prev(h2, mesh);
    halfedge_descriptor newh = halfedge(add_edge(mesh), mesh),
        newh_opp = opposite(newh, mesh);
    // set target vertices of the new halfedges
    set_target(newh, target(h1, mesh), mesh);
    set_target(newh_opp, target(ph2, mesh), mesh);
    // update next/prev pointers
    set_next(h1, newh_opp, mesh);
    set_next(newh_opp, h2, mesh);
    set_next(ph2, newh, mesh);
    set_next(newh, nh1, mesh);
  }
  for(std::size_t i = 0; i < input_halfedges.size(); ++i)
  {
    halfedge_descriptor h = input_halfedges[i];
    face_descriptor nf = add_face(mesh);

    std::array<halfedge_descriptor, 4> hedges;
    for (int k=0; k<4; ++k)
    {
      hedges[k]=h;
      h = next(h, mesh);
    }

    set_face(hedges[0], nf, mesh);
    set_face(hedges[1], nf, mesh);
    set_face(hedges[2], nf, mesh);
    set_face(hedges[3], nf, mesh);
    set_halfedge(nf, hedges[0], mesh);
    Euler::split_face(hedges[0], hedges[2], mesh);
  }
}

template<typename PMAP, typename Vector>
struct Const_dist_translation{
  Const_dist_translation(PMAP map, const Vector& dir)
    :map(map), dir(dir){}

  template<typename VertexDescriptor, typename U>
  void operator()(const VertexDescriptor vd, const U&) const
  {
    typename boost::property_traits<PMAP>::value_type p = get(map, vd) + dir;
    put(map, vd, p);
  }

  PMAP map;
  Vector dir;
};

struct Identity_functor
{
  template<typename T, typename U>
  void operator()(const T&, const U&) const {}
};
}//end extrude_impl

/**
 * \ingroup PMP_meshing_grp
 *
 * \brief performs a generalized extrusion of `input` and puts it in `output`.
 *
 * This function extrudes the open surface mesh `input` and puts the result in `output`. The mesh generated is a closed
 * surface mesh with a bottom and top part, both having the same graph combinatorics as `input` (except
 * that the orientation of the faces of the bottom part is reversed). The bottom and the top parts are
 * connected by a triangle strip between boundary cycles. The coordinates of the points associated to the
 * vertices of the bottom and top part are first initialized to the same value as the corresponding
 * vertices of `input`. Then for each vertex, a call to `bot` and `top` is done for the vertices of the
 * bottom part and the top part, respectively.
 *
 * \attention `output` may be self intersecting.
 *
 * @tparam InputMesh a model of `FaceListGraph`
 * @tparam OutputMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters" for `InputMesh`
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters" for `OutputMesh`
 * @tparam BottomFunctor a functor providing
 * \code {.cpp}
 * void operator()`(boost::graph_traits<InputMesh>::vertex_descriptor input_v,boost::graph_traits<OutputMesh>::vertex_descriptor output_v)
 * \endcode
 *  where `output_v` is the copy of `input_v` from `input` into the bottom part of `output`.
 *
 * @tparam TopFunctor a functor providing a similar `operator()` as `BottomFunctor`.
 *
 * @param input an open surface mesh to extrude.
 * @param output a surface mesh that will contain the result of the extrusion.
 * @param bot functor that will transform all points copied from
 * `input` in order to shape the bottom part of the extrusion.
 * @param top functor that will transform all points copied from
 * `input` in order to shape the top part of the extrusion.
 * @param np_in an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `input`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<InputMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, input)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `input`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `ouput`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<OutputMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, ouput)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `ouput`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template <class InputMesh,
          class OutputMesh,
          class BottomFunctor,
          class TopFunctor,
          class NamedParameters1,
          class NamedParameters2
          >
void extrude_mesh(const InputMesh& input,
                  OutputMesh& output,
                  const BottomFunctor& bot,
                  const TopFunctor& top,
                  const NamedParameters1& np_in,
                  const NamedParameters2& np_out)
{
  typedef typename boost::graph_traits<InputMesh>::vertex_descriptor input_vertex_descriptor;
  typedef typename boost::graph_traits<InputMesh>::halfedge_descriptor input_halfedge_descriptor;

  typedef typename boost::graph_traits<OutputMesh>::vertex_descriptor   output_vertex_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::halfedge_descriptor output_halfedge_descriptor;

  CGAL_assertion(!CGAL::is_closed(input));
  typedef typename GetVertexPointMap < OutputMesh, NamedParameters2>::type VPMap;
  typedef typename GetVertexPointMap < InputMesh, NamedParameters1>::const_type IVPMap;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  VPMap output_vpm = choose_parameter(get_parameter(np_out, internal_np::vertex_point),
                                      get_property_map(vertex_point, output));
  IVPMap input_vpm = choose_parameter(get_parameter(np_in, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, input));

  std::vector<std::pair<input_vertex_descriptor, output_vertex_descriptor> > bottom_v2v;
  std::vector<std::pair<input_halfedge_descriptor, output_halfedge_descriptor> > bottom_h2h;
  copy_face_graph(input, output, parameters::vertex_point_map(input_vpm)
                                            .vertex_to_vertex_output_iterator(std::back_inserter(bottom_v2v))
                                            .halfedge_to_halfedge_output_iterator(std::back_inserter(bottom_h2h)),
                                 parameters::vertex_point_map(output_vpm));

  // create the offset for the other side
  for(std::size_t i = 0; i< bottom_v2v.size(); ++i)
  {
    bot(bottom_v2v[i].first, bottom_v2v[i].second);
  }
  CGAL::Polygon_mesh_processing::reverse_face_orientations(output);

  // collect border halfedges for the creation of the triangle strip
  std::vector<std::pair<input_vertex_descriptor, output_vertex_descriptor> > top_v2v;
  std::vector<std::pair<input_halfedge_descriptor, output_halfedge_descriptor> > top_h2h;
  copy_face_graph(input, output, parameters::vertex_point_map(input_vpm)
                                            .vertex_to_vertex_output_iterator(std::inserter(top_v2v, top_v2v.end()))
                                            .halfedge_to_halfedge_output_iterator(std::inserter(top_h2h, top_h2h.end())),
                                 parameters::vertex_point_map(output_vpm));

  for(std::size_t i = 0; i< top_v2v.size(); ++i)
  {
    top(top_v2v[i].first, top_v2v[i].second);
  }
  std::vector<output_halfedge_descriptor> border_hedges;
  std::vector<output_halfedge_descriptor> offset_border_hedges;
  for(std::size_t i = 0; i< top_h2h.size(); ++i)
  {
    input_halfedge_descriptor h = top_h2h[i].first;
    if( CGAL::is_border(h, input) )
    {
      border_hedges.push_back(top_h2h[i].second);
      offset_border_hedges.push_back(bottom_h2h[i].second);
      CGAL_assertion(is_border(border_hedges.back(), output));
      CGAL_assertion(is_border(offset_border_hedges.back(), output));
    }
  }
  // now create a triangle strip
  extrude_impl::create_strip(border_hedges, offset_border_hedges, output);
}


/**
 * \ingroup PMP_meshing_grp
 *
 * fills `output` with a closed mesh bounding the volume swept by `input` when translating its
 * vertices by `v`. The mesh is oriented so that the faces corresponding to `input`
 * in `output` have the same orientation.
 *
 * \attention `output` may be self intersecting.
 *
 * @tparam InputMesh a model of the concept `FaceListGraph`
 * @tparam OutputMesh a model of the concept `FaceListGraph` and `MutableFaceGraph`
 * @tparam Vector_3 vector type from the same CGAL kernel as the point of the vertex point map used for `OutputMesh`.
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters" for `InputMesh`
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters" for `OutputMesh`
 *
 * @param input an open surface mesh to extrude.
 * @param output a surface mesh that will contain the result of the extrusion.
 * @param v the vector defining the direction of the extrusion
 * @param np_in an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `input`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<InputMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, input)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `input`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `output`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<OutputMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, output)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `output`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template <class InputMesh,
          class OutputMesh,
          class NamedParameters1,
          class NamedParameters2>
void extrude_mesh(const InputMesh& input,
                  OutputMesh& output,
                  #ifdef DOXYGEN_RUNNING
                  Vector_3 v,
                  #else
                  typename GetGeomTraits<OutputMesh, NamedParameters2>::type::Vector_3 v,
                  #endif
                  const NamedParameters1& np_in,
                  const NamedParameters2& np_out)
{
  typedef typename GetVertexPointMap < OutputMesh, NamedParameters2>::type VPMap;
  VPMap output_vpm = parameters::choose_parameter(parameters::get_parameter(np_out, internal_np::vertex_point),
                                  get_property_map(vertex_point, output));

  extrude_impl::Const_dist_translation<
      typename GetVertexPointMap<OutputMesh, NamedParameters2>::type,
      typename GetGeomTraits<OutputMesh, NamedParameters2>::type::Vector_3> bot(output_vpm,
                                                                                v);
  extrude_impl::Identity_functor top;
  extrude_mesh(input, output, bot,top, np_in, np_out);
}
//convenience overload
template <class InputMesh,
          class OutputMesh,
          typename Vector>
void extrude_mesh(const InputMesh& input,
                  OutputMesh& output,
                  Vector dir)
{
  extrude_mesh(input, output, dir,
               parameters::all_default(),
               parameters::all_default());
}

template <class InputMesh,
          class OutputMesh,
          typename Vector,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void extrude_mesh(const InputMesh& input,
                  OutputMesh& output,
                  Vector dir,
                  const CGAL_PMP_NP_CLASS& np)
{
  extrude_mesh(input, output, dir,
               np,
               parameters::all_default());
}

template <class InputMesh,
          class OutputMesh,
          class BottomFunctor,
          class TopFunctor>
void extrude_mesh(const InputMesh& input,
                  OutputMesh& output,
                  const BottomFunctor& bot,
                  const TopFunctor& top)
{
  extrude_mesh(input, output, bot, top,
               parameters::all_default(), parameters::all_default());
}

}} //end CGAL::PMP
#endif //CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H
