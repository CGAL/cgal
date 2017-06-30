// Copyright (c) 2015 GeometryFactory (France).
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
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
#define CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/unordered_set.hpp>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/algorithm.h>

#include <set>

namespace CGAL{
namespace Polygon_mesh_processing {

  namespace internal
  {
    template<typename PM
           , typename FaceRange
           , typename HalfedgeOutputIterator>
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& faces
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh)
    {
      typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
      typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

      //collect halfedges that appear only once
      // the bool is true if the halfedge stored is the one of the face,
      // false if it is its opposite
      std::map<halfedge_descriptor, bool> border;
      BOOST_FOREACH(face_descriptor f, faces)
      {
        BOOST_FOREACH(halfedge_descriptor h,
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
      BOOST_FOREACH(const HD_bool& hd, border)
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
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& faces
                                               , typename boost::cgal_no_property::type
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh
                                               , const NamedParameters& /* np */)
    {
      return border_halfedges_impl(faces, out, pmesh);
    }

    template<typename PM
           , typename FaceRange
           , typename FaceIndexMap
           , typename HalfedgeOutputIterator
           , typename NamedParameters>
    HalfedgeOutputIterator border_halfedges_impl(const FaceRange& faces
                                               , const FaceIndexMap& fmap
                                               , HalfedgeOutputIterator out
                                               , const PM& pmesh
                                               , const NamedParameters& /* np */)
    {
      typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
      typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;

      //make a minimal check that it's properly initialized :
      //if the 2 first faces have the same id, we know the property map is not initialized
      if (boost::is_same<typename GetFaceIndexMap<PM, NamedParameters>::Is_internal_map,
                         boost::true_type>::value)
      {
        typename boost::range_iterator<const FaceRange>::type it = boost::const_begin(faces);
        if (get(fmap, *it) == get(fmap, *cpp11::next(it)))
        {
          std::cerr << "WARNING : the internal property map for CGAL::face_index_t" << std::endl
                    << "          is not properly initialized." << std::endl
                    << "          Initialize it before calling border_halfedges()" << std::endl;
        }
      }

      std::vector<bool> present(num_faces(pmesh), false);
      BOOST_FOREACH(face_descriptor fd, faces)
        present[get(fmap, fd)] = true;

      BOOST_FOREACH(face_descriptor fd, faces)
        BOOST_FOREACH(halfedge_descriptor hd,
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
  \ingroup PkgPolygonMeshProcessing
  * collects the border halfedges of a surface patch defined as a face range.
  * For each returned halfedge `h`, `opposite(h, pmesh)` belongs to a face of the patch,
  * but `face(h, pmesh)` does not belong to the patch.
  *
  * @tparam PolygonMesh model of `HalfedgeGraph`. If `PolygonMesh`
  *  has an internal property map
  *  for `CGAL::face_index_t` and no `face_index_map` is given
  *  as a named parameter, then the internal one should be initialized
  * @tparam FaceRange range of
       `boost::graph_traits<PolygonMesh>::%face_descriptor`, model of `Range`.
        Its iterator type is `InputIterator`.
  * @tparam HalfedgeOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
     for patch border
  * @tparam NamedParameters a sequence of \ref namedparameters
  *
  * @param pmesh the polygon mesh to which `faces` belong
  * @param faces the range of faces defining the patch whose border halfedges
  *              are collected
  * @param out the output iterator that collects the border halfedges of the patch,
  *            seen from outside.
  * @param np optional sequence of \ref namedparameters among the ones listed below

  * \cgalNamedParamsBegin
      \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
    \cgalNamedParamsEnd
  *
  * @returns `out`
  */
  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator
         , typename NamedParameters>
  HalfedgeOutputIterator border_halfedges(const FaceRange& faces
                                  , const PolygonMesh& pmesh
                                  , HalfedgeOutputIterator out
                                  , const NamedParameters& np)
  {
    if (faces.empty()) return out;

    typedef PolygonMesh PM;
    typedef typename GetFaceIndexMap<PM, NamedParameters>::const_type     FIMap;
    typedef typename boost::property_map<typename internal::Dummy_PM,
                                              CGAL::face_index_t>::type   Unset_FIMap;

    if (boost::is_same<FIMap, Unset_FIMap>::value || faces.size() == 1)
    {
      //face index map is not given in named parameters, nor as an internal property map
      return internal::border_halfedges_impl(faces, out, pmesh);
    }

    //face index map given as a named parameter, or as an internal property map
    FIMap fim = boost::choose_param(get_param(np, internal_np::face_index),
                                    get_const_property_map(CGAL::face_index, pmesh));

    return internal::border_halfedges_impl(faces, fim, out, pmesh, np);
  }

  template<typename PolygonMesh
         , typename HalfedgeOutputIterator>
  HalfedgeOutputIterator border_halfedges(const PolygonMesh& pmesh
                                        , HalfedgeOutputIterator out)
  {
    typedef PolygonMesh PM;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(pmesh))
      if (is_border(hd, pmesh))
        *out++ = hd;
    return out;
  }

  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator>
  HalfedgeOutputIterator border_halfedges(const FaceRange& faces
                                        , const PolygonMesh& pmesh
                                        , HalfedgeOutputIterator out)
  {
    return border_halfedges(faces, pmesh, out,
      CGAL::Polygon_mesh_processing::parameters::all_default());
  }


  // counts the number of connected components of the boundary of the mesh.
  //
  // @tparam PolygonMesh model of `HalfedgeGraph`.
  //
  // @param pmesh the polygon mesh to which `faces` belong
  //
  template<typename PolygonMesh>
  unsigned int number_of_borders(const PolygonMesh& pmesh)
  {
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

    unsigned int border_counter = 0;
    boost::unordered_set<halfedge_descriptor> visited;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh)){
      if(visited.find(h)== visited.end()){
        if(is_border(h,pmesh)){
          ++border_counter;
          BOOST_FOREACH(halfedge_descriptor haf, halfedges_around_face(h, pmesh)){
            visited.insert(haf);
          }
        }
      }
    }

    return border_counter;
  }
} // end of namespace Polygon_mesh_processing
} // end of namespace CGAL


#endif //CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
