// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_BGL_BORDER_H
#define CGAL_BGL_BORDER_H

#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/graph/graph_traits.hpp>

#include <map>
#include <unordered_set>
#include <utility>

namespace CGAL {
namespace internal {

template <typename Graph>
std::size_t border_size(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                        const Graph& g)
{
  // if you want to use it on a non-border halfedge, just use `degree(face, mesh)`
  CGAL_precondition(is_border(h, g));

  std::size_t res = 0;
  typename boost::graph_traits<Graph>::halfedge_descriptor done = h;
  do
  {
    ++res;
    h = next(h, g);
  }
  while(h != done);

  return res;
}

template<typename Graph,
         typename FaceRange,
         typename HalfedgeOutputIterator>
HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range,
                                             HalfedgeOutputIterator out,
                                             const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor     face_descriptor;

  // collect halfedges that appear only once
  // the bool is true if the halfedge stored is the one of the face,
  // false if it is its opposite
  std::map<halfedge_descriptor, bool> border;
  for(face_descriptor f : face_range)
  {
    for(halfedge_descriptor h :
      halfedges_around_face(halfedge(f, g), g))
    {
      // halfedge_descriptor is model of `LessThanComparable`
      bool from_face = (h < opposite(h, g));
      halfedge_descriptor he = from_face ? h : opposite(h, g);
      if (border.find(he) != border.end())
        border.erase(he); //even number of appearances
      else
        border.insert(std::make_pair(he, from_face));//odd number of appearances
    }
  }
  //copy them in out
  typedef typename std::map<halfedge_descriptor, bool>::value_type HD_bool;
  for(const HD_bool& hd : border)
  {
    if (!hd.second) // to get the border halfedge (which is not on the face)
      *out++ = hd.first;
    else
      *out++ = opposite(hd.first, g);
  }
  return out;
}

template <typename Graph,
          typename FaceRange,
          typename HalfedgeOutputIterator,
          typename NamedParameters>
HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range,
                                             typename boost::cgal_no_property::type,
                                             HalfedgeOutputIterator out,
                                             const Graph& g,
                                             const NamedParameters& /* np */)
{
  return border_halfedges_impl(face_range, out, g);
}

template <typename Graph,
          typename FaceRange,
          typename FaceIndexMap,
          typename HalfedgeOutputIterator,
          typename NamedParameters>
HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range,
                                             const FaceIndexMap& fmap,
                                             HalfedgeOutputIterator out,
                                             const Graph& g,
                                             const NamedParameters& /* np */)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor     face_descriptor;

  CGAL_assertion(BGL::internal::is_index_map_valid(fmap, num_faces(g), faces(g)));

  std::vector<bool> present(num_faces(g), false);
  for(face_descriptor fd : face_range)
    present[get(fmap, fd)] = true;

  for(face_descriptor fd : face_range)
  {
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, g), g))
    {
      halfedge_descriptor opp = opposite(hd, g);
      if (is_border(opp, g) || !present[get(fmap, face(opp,g))])
        *out++ = opp;
    }
  }

  return out;
}

struct Dummy_PM
{
public:
  typedef bool vertex_property_type;
};

} // namespace internal

/*!
* \ingroup PkgBGLHelperFct
*
* \brief collects the border halfedges of a surface patch defined as a face range.
*
* For each returned halfedge `h`, `opposite(h, g)` belongs to a face of the patch,
* but `face(h, g)` does not belong to the patch.
*
* @tparam Graph model of `HalfedgeGraph`
* @tparam FaceRange a model of `Range` with value type `boost::graph_traits<Graph>::%face_descriptor`
* @tparam HalfedgeOutputIterator model of `OutputIterator` holding `boost::graph_traits<Graph>::%halfedge_descriptor` for patch border
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param g the polygon mesh to which the faces in `face_range` belong
* @param face_range the range of faces defining the patch whose border halfedges are collected
* @param out the output iterator that collects the border halfedges of the patch, seen from outside
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

* \cgalNamedParamsBegin
*   \cgalParamNBegin{face_index_map}
*     \cgalParamDescription{a property map associating to each face of `g` a unique index between `0` and `num_faces(g) - 1`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
*                    as key type and `std::size_t` as value type}
*     \cgalParamDefault{an automatically indexed internal map}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @returns `out`
*
* @see `extract_boundary_cycles()`
*/
template <typename Graph,
          typename FaceRange,
          typename HalfedgeOutputIterator,
          typename NamedParameters = parameters::Default_named_parameters>
HalfedgeOutputIterator border_halfedges(const FaceRange& face_range,
                                        const Graph& g,
                                        HalfedgeOutputIterator out,
                                        const NamedParameters& np = parameters::default_values())
{
  if (face_range.empty())
    return out;

  typedef typename CGAL::GetInitializedFaceIndexMap<Graph, NamedParameters>::const_type FIMap;
  FIMap fim = CGAL::get_initialized_face_index_map(g, np);

  return internal::border_halfedges_impl(face_range, fim, out, g, np);
}

template <typename Graph,
          typename HalfedgeOutputIterator>
HalfedgeOutputIterator border_halfedges(const Graph& g,
                                        HalfedgeOutputIterator out)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  for(halfedge_descriptor hd : halfedges(g))
    if(is_border(hd, g))
      *out++ = hd;
  return out;
}

// counts the number of connected components of the boundary of the mesh.
//
// @tparam Graph model of `HalfedgeGraph`.
//
// @param g the polygon mesh to which `face_range` belong
//
template<typename Graph>
unsigned int number_of_borders(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  unsigned int border_counter = 0;
  std::unordered_set<halfedge_descriptor> visited;
  for(halfedge_descriptor h : halfedges(g)){
    if(visited.find(h)== visited.end()){
      if(is_border(h,g)){
        ++border_counter;
        for(halfedge_descriptor haf : halfedges_around_face(h, g)){
          visited.insert(haf);
        }
      }
    }
  }

  return border_counter;
}

/*!
* @ingroup PkgBGLHelperFct
*
* extracts boundary cycles as a list of halfedges, with one halfedge per border.
*
* @tparam Graph a model of `HalfedgeListGraph`
* @tparam OutputIterator a model of `OutputIterator` holding objects of type
*   `boost::graph_traits<Graph>::%halfedge_descriptor`
*
* @param g a polygon mesh
* @param out an output iterator where the border halfedges will be put
*
* @see `border_halfedges()`
*/
template <typename Graph, typename OutputIterator>
OutputIterator extract_boundary_cycles(const Graph& g,
                                        OutputIterator out)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  std::unordered_set<halfedge_descriptor> hedge_handled;
  for(halfedge_descriptor h : halfedges(g))
  {
    if(is_border(h, g) && hedge_handled.insert(h).second)
    {
      *out++ = h;
      for(halfedge_descriptor h2 : halfedges_around_face(h, g))
        hedge_handled.insert(h2);
    }
  }
  return out;
}

} // namespace CGAL

#endif // CGAL_BGL_BORDER_H
