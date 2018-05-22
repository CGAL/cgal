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


#ifndef CGAL_EXTUR_POLYGON_MESH_H
#define CGAL_EXTUR_POLYGON_MESH_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/Kernel_traits.h>
#include <boost/unordered_map.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {

/**
 * \ingroup ???
 * Extrudes `input` into `output` following the direction given by `dir` and
 * at a distance `d`.
 * @tparam InputMesh a model of the concept `FaceListGraph`
 * @tparam OutputMesh a model of the concept `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 * @param input the closed triangulated `InputMesh` to extrude.
 * @param output the `OutputMesh` containing the result of the extrusion.
 * @param np an optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
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
          class NamedParameters>
void extrude_mesh(const InputMesh& input, 
                  OutputMesh& output, 
                  #ifdef DOXYGEN_RUNNING
                  Vector_3 dir,
                  #else
                  typename CGAL::Kernel_traits<
                  typename boost::property_traits<
                  typename GetVertexPointMap <
                  OutputMesh, NamedParameters>::type>::
                  value_type>::Kernel::Vector_3 dir, 
                  #endif
                  const double d, 
                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<InputMesh>::halfedge_descriptor input_halfedge_descriptor;
  
  typedef typename boost::graph_traits<OutputMesh>::vertex_descriptor   output_vertex_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::halfedge_descriptor output_halfedge_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::edge_descriptor     output_edge_descriptor;
  typedef typename boost::graph_traits<OutputMesh>::face_descriptor     output_face_descriptor;
  
  CGAL_assertion(!CGAL::is_closed(input));
  typedef typename GetVertexPointMap < OutputMesh, NamedParameters>::type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type Point;
  
  boost::unordered_map<input_halfedge_descriptor, output_halfedge_descriptor> offset_h2h;
  copy_face_graph(input, output, CGAL::Emptyset_iterator(),
                        std::inserter(offset_h2h, offset_h2h.end()));
  VPMap output_vpm = choose_param(get_param(np, internal_np::vertex_point),
                                   get_property_map(vertex_point, output));
  // create the offset for the other side
  BOOST_FOREACH(output_vertex_descriptor v, vertices(output))
  {
    Point p = get(output_vpm, v) + d*dir;
    put(output_vpm, v, p);
  }
  CGAL::Polygon_mesh_processing::reverse_face_orientations(output);
  
  // collect border halfedges for the creation of the triangle strip
  boost::unordered_map<input_halfedge_descriptor, output_halfedge_descriptor> h2h;
  copy_face_graph(input, output, CGAL::Emptyset_iterator(),
                        std::inserter(h2h, h2h.end()));
  std::vector<output_halfedge_descriptor> border_hedges;
  std::vector<output_halfedge_descriptor> offset_border_hedges;
  BOOST_FOREACH(input_halfedge_descriptor h, halfedges(input))
  {
    if( CGAL::is_border(h, input) )
    {
      border_hedges.push_back(h2h[h]);
      offset_border_hedges.push_back(offset_h2h[h]);
      CGAL_assertion(is_border(border_hedges.back(), output));
      CGAL_assertion(is_border(offset_border_hedges.back(), output));
    }
  }
  // now create a quad strip
  for(std::size_t i=0; i< border_hedges.size(); ++i)
  {
    //     before                 after
    // -----  o  -------     -----  o  -------
    // <----     <-----      <----  |   <-----
    //  nh1        h1         nh1   |     h1
    //                              |
    //                        hnew  |  hnew_opp
    //                              |
    //   ph2       h2          ph2  |     h2
    //  ---->    ----->       ----> |   ----->
    // -----  o  -------     -----  o  -------
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
    
    output_face_descriptor nf1 = add_face(output);
    output_face_descriptor nf2 = add_face(output);
    
    CGAL::cpp11::array<output_halfedge_descriptor, 4> hedges;
    for (int k=0; k<4; ++k)
    {
      hedges[k]=h;
      h = next(h, output);
    }
    
    //add a diagonal
    output_edge_descriptor new_e = add_edge(output);
    output_halfedge_descriptor new_h = halfedge(new_e, output),
        new_h_opp = opposite(new_h, output);
    // set vertex pointers
    set_target(new_h_opp, target(hedges[0], output), output);
    set_target(new_h, target(hedges[2], output), output);
    
    // set next pointers
    set_next(hedges[0], new_h, output);
    set_next(new_h, hedges[3], output);
    set_next(hedges[2], new_h_opp, output);
    set_next(new_h_opp, hedges[1], output);
    
    // set face halfedge pointers
    set_face(hedges[0], nf1, output);
    set_face(hedges[3], nf1, output);
    set_face(new_h, nf1, output);
    set_face(hedges[1], nf2, output);
    set_face(hedges[2], nf2, output);
    set_face(new_h_opp, nf2, output);
    
    // set halfedge face pointers
    set_halfedge(nf1, hedges[0], output);
    set_halfedge(nf2, hedges[2], output);
  }
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
               parameters::all_default());
}

}} //end CGAL::PMP
#endif //CGAL_EXTRUDE_POLYGON_MESH_H
