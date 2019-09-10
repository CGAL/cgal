// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H
#define CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/unordered.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/array.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationMesh

    \brief Class that  generates models of `NeighborQuery` based on
    an input mesh.

    \tparam FaceListGraph model of `FaceListGraph`. 
  */
template <typename FaceListGraph>
class Mesh_neighborhood
{
public:
  typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor; ///<

private:
  typedef typename boost::graph_traits<FaceListGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceListGraph>::vertex_descriptor vertex_descriptor;
  const FaceListGraph& m_mesh;

  class Is_face_selected
  {
  public:
    typedef face_descriptor key_type;
    typedef bool value_type;
    typedef bool reference;
    typedef boost::read_write_property_map_tag category;
    typedef typename CGAL::cpp11::unordered_set<face_descriptor, CGAL::Handle_hash_function> Set;
  private:
    Set* m_set;
    
  public:
    Is_face_selected(Set* set = NULL) : m_set (set) { }

    inline friend value_type get (const Is_face_selected& pm, const key_type& k)
    {
      return (pm.m_set->find(k) != pm.m_set->end());
    }
   
    inline friend void put (const Is_face_selected& pm, const key_type& k, const value_type&)
    {
      pm.m_set->insert(k);
    }
  };
  
public:

  /*!
    Functor that computes the 1-ring neighborhood of the face of an input mesh.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Mesh_neighborhood
  */
  class One_ring_neighbor_query
  {
  public:
    typedef typename Mesh_neighborhood::face_descriptor value_type; ///<
  private:
    const Mesh_neighborhood& neighborhood;

  public:

    /*!
      \brief Constructs a 1-ring neighbor query object.
      \param neighborhood mesh neighborhood object.
    */
    One_ring_neighbor_query (const Mesh_neighborhood& neighborhood)
      : neighborhood (neighborhood) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.one_ring_neighbors (query, output);
      return output;
    }
    /// \endcond
  };

  /*!
    Functor that computes the N-ring neighborhood of the face of an input mesh.

    \cgalModels CGAL::Classification::NeighborQuery

    \sa Mesh_neighborhood
  */
  class N_ring_neighbor_query
  {
  public:
    typedef typename Mesh_neighborhood::face_descriptor value_type; ///<
  private:
    const Mesh_neighborhood& neighborhood;
    const std::size_t n;

  public:

    /*!
      \brief Constructs a N-ring neighbor query object.
      \param neighborhood mesh neighborhood object.
      \param n size of neighborhood.
    */
    N_ring_neighbor_query (const Mesh_neighborhood& neighborhood, const std::size_t n)
      : neighborhood (neighborhood), n(n) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.n_ring_neighbors (query, output, n);
      return output;
    }
    /// \endcond
  };

  /// \cond SKIP_IN_MANUAL
  friend class One_ring_neighbor_query;
  friend class N_ring_neighbor_query;
  /// \endcond

  /// \name Constructor
  /// @{

  /*!
    \brief Constructs a neighborhood object based on the input mesh.

    \param mesh input mesh.
  */
  Mesh_neighborhood (const FaceListGraph& mesh) : m_mesh (mesh)
  {
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  ~Mesh_neighborhood ()
  {
  }
  /// \endcond

  /// \name Queries
  /// @{

  /*!
    \brief Returns a 1-ring neighbor query object.
  */
  One_ring_neighbor_query one_ring_neighbor_query () const
  {
    return One_ring_neighbor_query (*this);
  }

  /*!
    \brief Returns an N-ring neighbor query object.
  */
  N_ring_neighbor_query n_ring_neighbor_query (const std::size_t n) const
  {
    return N_ring_neighbor_query (*this, n);
  }

  /// @}


private:

  template <typename OutputIterator>
  void direct_neighbors (const face_descriptor& query, OutputIterator output) const
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(query, m_mesh), m_mesh))
      {
        *(output ++ ) = face(opposite(hd, m_mesh), m_mesh);
      }
  }
  
  template <typename OutputIterator>
  void one_ring_neighbors (const face_descriptor& query, OutputIterator output) const
  {
    return n_ring_neighbors (query, output, 1);
  }

  template <typename OutputIterator>
  void n_ring_neighbors (const face_descriptor& query, OutputIterator output, const std::size_t n) const
  {
    *(output ++) = get(get(CGAL::face_index, m_mesh), query);
    CGAL::cpp11::array<face_descriptor,1> init = {{ query }};
    typename Is_face_selected::Set done;
    done.insert(query);
    std::vector<face_descriptor> desc;
    CGAL::expand_face_selection
      (init, m_mesh, static_cast<unsigned int>(n), Is_face_selected(&done), std::back_inserter (desc));
    for (std::size_t i = 0; i < desc.size(); ++ i)
      *(output ++) = get(get(CGAL::face_index, m_mesh), desc[i]);
  }


};
  

}
  
}


#endif // CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H
