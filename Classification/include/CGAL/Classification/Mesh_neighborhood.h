// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H
#define CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H

#include <CGAL/license/Classification.h>

#include <vector>

#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/properties.h>

#include <unordered_set>

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
  using face_descriptor = typename boost::graph_traits<FaceListGraph>::face_descriptor; ///<

private:
  using halfedge_descriptor = typename boost::graph_traits<FaceListGraph>::halfedge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<FaceListGraph>::vertex_descriptor;
  const FaceListGraph& m_mesh;

  class Is_face_selected
  {
  public:
    using key_type = face_descriptor;
    using value_type = bool;
    using reference = bool;
    using category = boost::read_write_property_map_tag;
    using Set = typename std::unordered_set<face_descriptor, CGAL::Handle_hash_function>;
  private:
    Set* m_set;

  public:
    Is_face_selected(Set* set = nullptr) : m_set (set) { }

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
    using value_type = typename Mesh_neighborhood::face_descriptor; ///<
  private:
    const Mesh_neighborhood& neighborhood;

  public:

    /*!
      \brief constructs a 1-ring neighbor query object.
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
    using value_type = typename Mesh_neighborhood::face_descriptor; ///<
  private:
    const Mesh_neighborhood& neighborhood;
    const std::size_t n;

  public:

    /*!
      \brief constructs a N-ring neighbor query object.
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
    \brief constructs a neighborhood object based on the input mesh.

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
    \brief returns a 1-ring neighbor query object.
  */
  One_ring_neighbor_query one_ring_neighbor_query () const
  {
    return One_ring_neighbor_query (*this);
  }

  /*!
    \brief returns an N-ring neighbor query object.
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
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(query, m_mesh), m_mesh))
      *(output ++ ) = face(opposite(hd, m_mesh), m_mesh);
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
    std::array<face_descriptor,1> init = {{ query }};
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
