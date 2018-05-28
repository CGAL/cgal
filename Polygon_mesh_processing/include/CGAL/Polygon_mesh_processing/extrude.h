// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot, Maxime Gimeno


#ifndef CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H
#define CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace extrude_impl{

template<typename PMAP, typename Vector>
struct Const_dist_translation{
  Const_dist_translation(PMAP map, const Vector& dir, const double d)
    :map(map), dir(dir), d(d){}
  
  template<typename VertexDescriptor, typename U>
  void operator()(const VertexDescriptor vd, const U&)
  {
    typename boost::property_traits<PMAP>::value_type p = get(map, vd) + d*dir;
    put(map, vd, p);
  }
  
  PMAP map;
  Vector dir;
  double d;
};

struct Identity_functor
{
  template<typename T, typename U>
  void operator()(const T&, const U&){}
};
}//end extrude_impl

/**
 * \ingroup PMP_meshing_grp
 * extrudes the open surface mesh `input` and puts the result in `output`. The mesh generated is a closed 
 * surface mesh with a bottom and top part, both having the same graph combinatorics as `input` (except 
 * that the orientation of the faces of the bottom part is reversed). The bottom and the top parts are 
 * connected by a triangle strip between each boundary cycles. The coordinates of the points associated to the 
 * vertices of the bottom and top part are first initialized to the same value as the corresponding 
 * vertices of `input`. Then for each vertex, a call to `bot` and `top` is done for the vertices of the 
 * bottom part and the top part, respectively.
 * @tparam InputMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * @tparam OutputMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters" for `InputMesh`
 * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters" for `OutputMesh`
 * @tparam BottomFunctor a functor with a function 
 * `void operator()(boost::graph_traits<InputMesh>::vertex_descriptor input_v,boost::graph_traits<OutputMesh>::vertex_descriptor output_v)`
 * , where `output_v` is the copy of `input_v` into the bottom part.
 * 
 * @tparam TopFunctor a functor with a an operator() similar to `BottomFunctor`.
 * @param input the open `InputMesh` to extrude.
 * @param output the `OutputMesh` containing the result of the extrusion.
 * @param bot a `BottomFunctor` that will apply a transformation to all points of 
 * `input` in order to create the first offsetted part of the extrusion. 
 * @param top a `TopFunctor` that will apply a transformation to all points of 
 * `input` in order to create the second offsetted part of the extrusion. 
 * @param np_in an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map that contains the points associated to the vertices of `input`. 
 * If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` 
 * should be available for the vertices of `input` \cgalParamEnd
 * \cgalNamedParamsEnd
 * 
 * * @param np_out an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map}
 *    the property map that will contain the points associated to the vertices of `output`. 
 * If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` 
 * should be available for the vertices of `output` \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <class InputMesh,
          class OutputMesh,
          class BottomFunctor,
          class TopFunctor,
          class NamedParameters1,
          class NamedParameters2
          >
void generic_extrude_mesh(const InputMesh& input, 
                          OutputMesh& output, 
                          BottomFunctor& bot,
                          TopFunctor& top,
                          const NamedParameters1& np_in,
                          const NamedParameters2& np_out)
{
  typedef typename boost::graph_traits<InputMesh>::vertex_descriptor input_vertex_descriptor;
  typedef typename boost::graph_traits<InputMesh>::halfedge_descriptor input_halfedge_descriptor;
  
  typedef typename boost::graph_traits<OutputMesh>::vertex_descriptor   output_vertex_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::halfedge_descriptor output_halfedge_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::face_descriptor     output_face_descriptor;
  
  CGAL_assertion(!CGAL::is_closed(input));
  typedef typename GetVertexPointMap < OutputMesh, NamedParameters2>::type VPMap;
  typedef typename GetVertexPointMap < InputMesh, NamedParameters1>::const_type IVPMap;
  
  VPMap output_vpm = choose_param(get_param(np_out, internal_np::vertex_point),
                                   get_property_map(vertex_point, output));
  IVPMap input_vpm = choose_param(get_param(np_in, internal_np::vertex_point),
                                   get_const_property_map(vertex_point, input));
  
  std::vector<std::pair<input_vertex_descriptor, output_vertex_descriptor> > bottom_v2v;
  std::vector<std::pair<input_halfedge_descriptor, output_halfedge_descriptor> > bottom_h2h;
  copy_face_graph(input, output, std::back_inserter(bottom_v2v),
                        std::back_inserter(bottom_h2h), Emptyset_iterator(),
                  input_vpm, output_vpm);
  
  // create the offset for the other side
  for(std::size_t i = 0; i< bottom_v2v.size(); ++i)
  {
    bot(bottom_v2v[i].first, bottom_v2v[i].second);
  }
  CGAL::Polygon_mesh_processing::reverse_face_orientations(output);
  
  // collect border halfedges for the creation of the triangle strip
  std::vector<std::pair<input_vertex_descriptor, output_vertex_descriptor> > top_v2v;
  std::vector<std::pair<input_halfedge_descriptor, output_halfedge_descriptor> > top_h2h;
  copy_face_graph(input, output, std::inserter(top_v2v, top_v2v.end()),
                        std::inserter(top_h2h, top_h2h.end()), Emptyset_iterator(),
                  input_vpm, output_vpm);
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
  for(std::size_t i=0; i< border_hedges.size(); ++i)
  {
    output_halfedge_descriptor h1 = border_hedges[i], h2=offset_border_hedges[i],
        nh1 = next(h1, output), ph2 = prev(h2, output);
    output_halfedge_descriptor newh = halfedge(add_edge(output), output),
        newh_opp = opposite(newh, output);
    // set target vertices of the new halfedges
    set_target(newh, target(h1, output), output);
    set_target(newh_opp, target(ph2, output), output);
    // update next/prev pointers
    set_next(h1, newh_opp, output);
    set_next(newh_opp, h2, output);
    set_next(ph2, newh, output);
    set_next(newh, nh1, output);
  }
  
  // create new faces
  for(std::size_t i=0; i< border_hedges.size(); ++i)
  {
    output_halfedge_descriptor h = border_hedges[i];
    
    output_face_descriptor nf = add_face(output);
    
    CGAL::cpp11::array<output_halfedge_descriptor, 4> hedges;
    for (int k=0; k<4; ++k)
    {
      hedges[k]=h;
      h = next(h, output);
    }
    
    set_face(hedges[0], nf, output);
    set_face(hedges[1], nf, output);
    set_face(hedges[2], nf, output);
    set_face(hedges[3], nf, output);
    set_halfedge(nf, hedges[0], output);
    Euler::split_face(hedges[0], hedges[2], output);
  }
}


/**
 * \ingroup PMP_meshing_grp
 * fills `output` with a close mesh bounding the volume swept by `input` when translating its 
 * vertices by `dir` * `d`. The mesh is oriented so that the faces corresponding to `input`
 * in `output` have the same orientation.
 * @tparam InputMesh a model of the concept `MutableFaceGraph`
 * @tparam OutputMesh a model of the concept `MutableFaceGraph`
 * @tparam Vector_3 a type of Vector_3 from the kernel used by `OutputMesh`.
 * @tparam FT a type of floating type from the kernel used by `OutputMesh`.
 * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters" for `InputMesh`
 * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters" for `OutputMesh`
 * @tparam InputMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * @tparam OutputMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * @param input the open `InputMesh` to extrude.
 * @param output the `OutputMesh` containing the result of the extrusion.
 * @param dir the vector defining the direction of the extrusion
 * @param d the distance of the extrusion
 * @param np_in an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map}
 * the property map that contains the points associated to the vertices of `input`. 
 * If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` 
 * should be available for the vertices of `input` \cgalParamEnd
 * \cgalNamedParamsEnd
 * 
 * * @param np_out an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 * 
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map}
 * the property map that will contain the points associated to the vertices of `output`. 
 * If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` 
 * should be available for the vertices of `output` \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <class InputMesh,
          class OutputMesh,
          class NamedParameters1,
          class NamedParameters2>
void extrude_mesh(const InputMesh& input, 
                  OutputMesh& output, 
                  #ifdef DOXYGEN_RUNNING
                  Vector_3 dir,
                  const FT d, 
                  #else
                  typename GetGeomTraits<OutputMesh, NamedParameters2>::type::Vector_3 dir, 
                  const typename GetGeomTraits<OutputMesh, NamedParameters2>::type::FT d,
                  #endif
                  const NamedParameters1& np_in,
                  const NamedParameters2& np_out)
{
  typedef typename GetVertexPointMap < OutputMesh, NamedParameters2>::type VPMap;
  VPMap output_vpm = choose_param(get_param(np_out, internal_np::vertex_point),
                                   get_property_map(vertex_point, output));
  
  extrude_impl::Const_dist_translation<
      typename GetVertexPointMap<OutputMesh, NamedParameters2>::type,
      typename GetGeomTraits<OutputMesh, NamedParameters2>::type::Vector_3> bot(output_vpm, 
                                                                                  dir, d);
  extrude_impl::Identity_functor top;
  generic_extrude_mesh(input, output, bot,top, np_in, np_out);
}
//convenience overload
template <class InputMesh,
          class OutputMesh,
          typename Vector>
void extrude_mesh(const InputMesh& input, 
                  OutputMesh& output, 
                  Vector dir, 
                  const double d)
{
  extrude_mesh(input, output, dir, d, 
               parameters::all_default(),
               parameters::all_default());
}

template <class InputMesh,
          class OutputMesh,
          class BottomFunctor,
          class TopFunctor>
void generic_extrude_mesh(const InputMesh& input, 
                          OutputMesh& output, 
                          BottomFunctor& bot,
                          TopFunctor& top)
{
  generic_extrude_mesh(input, output, bot, top, 
                       parameters::all_default(), parameters::all_default());
}

template <class InputMesh,
          class OutputMesh,
          class BottomFunctor>
void generic_extrude_mesh(const InputMesh& input, 
                          OutputMesh& output, 
                          BottomFunctor& bot)
{
  extrude_impl::Identity_functor top;
  generic_extrude_mesh(input, output, bot, top, 
                       parameters::all_default(), parameters::all_default());
}

}} //end CGAL::PMP
#endif //CGAL_POLYGON_MESH_PROCESSING_EXTRUDE_H
