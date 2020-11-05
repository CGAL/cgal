// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
#define CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/unordered_set.hpp>

#include <set>

namespace CGAL{
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh>
std::size_t border_size(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                        const PolygonMesh& pmesh)
{
  // if you want to use it on a non-border halfedge, just use `degree(face, mesh)`
  CGAL_precondition(is_border(h, pmesh));

  std::size_t res = 0;
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor done = h;
  do
  {
    ++res;
    h = next(h, pmesh);
  }
  while(h != done);

  return res;
}

    template<typename PM
           , typename FaceRange
           , typename HalfedgeOutputIterator>
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh)
    {
      typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
      typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

      //collect halfedges that appear only once
      // the bool is true if the halfedge stored is the one of the face,
      // false if it is its opposite
      std::map<halfedge_descriptor, bool> border;
      for(face_descriptor f : face_range)
      {
        for(halfedge_descriptor h :
          halfedges_around_face(halfedge(f, pmesh), pmesh))
        {
          //halfedge_descriptor is model of `LessThanComparable`
          bool from_face = (h < opposite(h, pmesh));
          halfedge_descriptor he = from_face ? h : opposite(h, pmesh);
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
          *out++ = opposite(hd.first, pmesh);
      }
      return out;
    }

    template<typename PM
           , typename FaceRange
           , typename HalfedgeOutputIterator
           , typename NamedParameters>
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range
                                               , typename boost::cgal_no_property::type
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh
                                               , const NamedParameters& /* np */)
    {
      return border_halfedges_impl(face_range, out, pmesh);
    }

    template<typename PM
           , typename FaceRange
           , typename FaceIndexMap
           , typename HalfedgeOutputIterator
           , typename NamedParameters>
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& face_range
                                               , const FaceIndexMap& fmap
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh
                                               , const NamedParameters& /* np */)
    {
      typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
      typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

      CGAL_assertion(BGL::internal::is_index_map_valid(fmap, num_faces(pmesh), faces(pmesh)));

      std::vector<bool> present(num_faces(pmesh), false);
      for(face_descriptor fd : face_range)
        present[get(fmap, fd)] = true;

      for(face_descriptor fd : face_range)
        for(halfedge_descriptor hd :
                      halfedges_around_face(halfedge(fd, pmesh), pmesh))
       {
         halfedge_descriptor opp=opposite(hd, pmesh);
         if (is_border(opp, pmesh) || !present[get(fmap,face(opp,pmesh))])
           *out++ = opp;
       }

      return out;
    }

    struct Dummy_PM
    {
    public:
      typedef bool vertex_property_type;
    };

  }//end namespace internal

  /*!
  \ingroup PkgPolygonMeshProcessingRef
  * collects the border halfedges of a surface patch defined as a face range.
  * For each returned halfedge `h`, `opposite(h, pmesh)` belongs to a face of the patch,
  * but `face(h, pmesh)` does not belong to the patch.
  *
  * @tparam PolygonMesh model of `HalfedgeGraph`
  * @tparam FaceRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%face_descriptor`.
  * @tparam HalfedgeOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
     for patch border
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh to which the faces in `face_range` belong
  * @param face_range the range of faces defining the patch whose border halfedges
  *                   are collected
  * @param out the output iterator that collects the border halfedges of the patch,
  *            seen from outside.
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
  *                    as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @returns `out`
  */
  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator
         , typename NamedParameters>
  HalfedgeOutputIterator border_halfedges(const FaceRange& face_range
                                  , const PolygonMesh& pmesh
                                  , HalfedgeOutputIterator out
                                  , const NamedParameters& np)
  {
    if (face_range.empty())
      return out;

    typedef typename CGAL::GetInitializedFaceIndexMap<PolygonMesh, NamedParameters>::const_type FIMap;
    FIMap fim = CGAL::get_initialized_face_index_map(pmesh, np);

    return internal::border_halfedges_impl(face_range, fim, out, pmesh, np);
  }

  template<typename PolygonMesh
         , typename HalfedgeOutputIterator>
  HalfedgeOutputIterator border_halfedges(const PolygonMesh& pmesh
                                        , HalfedgeOutputIterator out)
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    for(halfedge_descriptor hd : halfedges(pmesh))
      if (is_border(hd, pmesh))
        *out++ = hd;
    return out;
  }

  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator>
  HalfedgeOutputIterator border_halfedges(const FaceRange& face_range
                                        , const PolygonMesh& pmesh
                                        , HalfedgeOutputIterator out)
  {
    return border_halfedges(face_range, pmesh, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }


  // counts the number of connected components of the boundary of the mesh.
  //
  // @tparam PolygonMesh model of `HalfedgeGraph`.
  //
  // @param pmesh the polygon mesh to which `face_range` belong
  //
  template<typename PolygonMesh>
  unsigned int number_of_borders(const PolygonMesh& pmesh)
  {
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    unsigned int border_counter = 0;
    boost::unordered_set<halfedge_descriptor> visited;
    for(halfedge_descriptor h : halfedges(pmesh)){
      if(visited.find(h)== visited.end()){
        if(is_border(h,pmesh)){
          ++border_counter;
          for(halfedge_descriptor haf : halfedges_around_face(h, pmesh)){
            visited.insert(haf);
          }
        }
      }
    }

    return border_counter;
  }

  /// @ingroup PkgPolygonMeshProcessingRef
  /// extracts boundary cycles as a list of halfedges, with one halfedge per border.
  ///
  /// @tparam PolygonMesh a model of `HalfedgeListGraph`
  /// @tparam OutputIterator a model of `OutputIterator` holding objects of type
  ///   `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
  ///
  /// @param pm a polygon mesh
  /// @param out an output iterator where the border halfedges will be put
  ///
  /// @todo It could make sense to also return the length of each cycle.
  /// @todo It should probably go into BGL package (like the rest of this file).
  template <typename PolygonMesh, typename OutputIterator>
  OutputIterator extract_boundary_cycles(PolygonMesh& pm,
                                         OutputIterator out)
  {
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    boost::unordered_set<halfedge_descriptor> hedge_handled;
    for(halfedge_descriptor h : halfedges(pm))
    {
      if(is_border(h, pm) && hedge_handled.insert(h).second)
      {
        *out++ = h;
        for(halfedge_descriptor h2 : halfedges_around_face(h, pm))
          hedge_handled.insert(h2);
      }
    }
    return out;
  }

} // end of namespace Polygon_mesh_processing
} // end of namespace CGAL


#endif //CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
