// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_TRIANGULATION_CONFORMER_2_H
#define CGAL_TRIANGULATION_CONFORMER_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_2/Refine_edges_with_clusters.h>

namespace CGAL {

template <typename Tr>
class Triangulation_conformer_2
{
  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;

  typedef Mesh_2::Refine_edges_with_clusters<Tr,
    Mesh_2::Is_locally_conforming_Gabriel<Tr> > Edges_level_Gabriel;

  typedef Mesh_2::Refine_edges_with_clusters<Tr,
    Mesh_2::Is_locally_conforming_Delaunay<Tr> > Edges_level_Delaunay;

protected:
  /** \name INITIALIZED */

  enum Initialization {
    NONE,     /**< \c this is not initialized. */
    CLUSTERS, /**< \c this clusters are initialized. */
    DELAUNAY, /**< \c this has been \e Delaunay-initialized. */
    GABRIEL   /**< \c this has been \e Gabriel-initialized. */
  };

// --- PROTECTED DATAS ---
  Initialization initialized;
  Tr& tr;
  Null_mesher_level null_level;
  Null_mesh_visitor null_visitor;
  Mesh_2::Clusters<Tr> clusters;
  Edges_level_Gabriel edges_level_Gabriel;
  Edges_level_Delaunay edges_level_Delaunay;

public:
  /** \name CONSTRUCTORS */
  Triangulation_conformer_2(Tr& tr_)
    : initialized(NONE),
      tr(tr_),
      null_level(), null_visitor(),
      clusters(tr_),
      edges_level_Gabriel(tr, clusters, null_level),
      edges_level_Delaunay(tr, clusters, null_level)
  {
  }

private:
  /** \name CHECKING METHODS*/

  template <typename Is_locally_conforming>
  bool is_conforming_XXX(Is_locally_conforming is_locally_conforming) const
  {
    for(Finite_edges_iterator ei = tr.finite_edges_begin();
        ei != tr.finite_edges_end();
        ++ei)
      if(ei->first->is_constrained(ei->second) &&
         !is_locally_conforming(tr, ei->first, ei->second) )
        return false;
    return true;
  }

public:  /** \name ACCESS TO CLUSTERS */
  typedef typename Mesh_2::Clusters<Tr>::Cluster_vertices_iterator
    Cluster_vertices_iterator;
  typedef typename Mesh_2::Clusters<Tr>::Vertices_in_cluster_iterator
    Vertices_in_cluster_iterator;

public:
  /** \name ACCESS FUNCTIONS */

  /** Access to the private initialized member data. */
  //@{
  void set_initialized(Initialization init) { initialized = init; }
  Initialization get_initialized() const { return initialized; }
  //@}

  int number_of_constrained_edges()
  {
    int nedges = 0;
    for(typename Tr::Finite_edges_iterator eit = tr.finite_edges_begin();
        eit != tr.finite_edges_end();
        ++eit)
      if(eit->first->is_constrained(eit->second))
        ++nedges;
    return nedges;
  }

  int number_of_clusters_vertices() const
  {
    return clusters.size();
  }

  Cluster_vertices_iterator clusters_vertices_begin() const
  {
    return clusters.clusters_vertices_begin();
  }

  Cluster_vertices_iterator clusters_vertices_end() const
  {
    return clusters.clusters_vertices_end();
  }

  unsigned int number_of_clusters_at_vertex(const Vertex_handle& vh)
  {
    return clusters.number_of_clusters_at_vertex(vh);
  }

#if 0
  // returns the sequence of vertices belonging to the n-th cluster of vh
  std::pair<Vertices_in_cluster_iterator, Vertices_in_cluster_iterator>
  vertices_in_cluster_sequence(const Vertex_handle& vh,
                               const unsigned int n)
  {
    return clusters.vertices_in_cluster_sequence();
  }
#endif

public:
  /** \name CHECKING METHODS */

  bool is_conforming_Delaunay()
  {
    typedef typename Mesh_2::Is_locally_conforming_Delaunay<Tr> Is_loc_conf;

    return is_conforming_XXX(Is_loc_conf());
  }

  bool is_conforming_Gabriel()
  {
    typedef typename Mesh_2::Is_locally_conforming_Gabriel<Tr> Is_loc_conf;

    return is_conforming_XXX(Is_loc_conf());
  }

  /** \name CONFORMING FUNCTIONS */

  void make_conforming_Delaunay()
  {
    if(initialized!=DELAUNAY) init_Delaunay();
    edges_level_Delaunay.refine(null_visitor);
  }

  void make_conforming_Gabriel()
  {
    if(initialized!=GABRIEL) init_Gabriel();
    edges_level_Gabriel.refine(null_visitor);
  }

  /** \name STEP BY STEP FUNCTIONS */

  // Note: step by step functions are not efficient at all!
private:
  void init_clusters()
  {
    if(initialized == NONE)
      clusters.create_clusters();
    initialized = CLUSTERS;
  }

public:
  /**
     Initializes the data structures
     (The call of this function is REQUIRED before any step by step
     operation).
  */
  //@{
  void init_Delaunay()
    {
      init_clusters();
      initialized = DELAUNAY;
      edges_level_Delaunay.scan_triangulation();
    }
  void init_Gabriel()
    {
      init_clusters();
      initialized = GABRIEL;
      edges_level_Gabriel.scan_triangulation();
    }
  //@}

  /** Tells if all constrained edges are conformed. */
  bool is_conforming_done()
    // This function cannot be "const" because, as edges_to_be_conformed is
    // filtred, its empty() method is not const.
  { return ( edges_level_Gabriel.no_longer_element_to_refine()
             && edges_level_Delaunay.no_longer_element_to_refine() );
  }

  /** Execute on step of the algorithm.
      init_XXX() should have been called before.
  */
  //@{
  bool try_one_step_conforming_Delaunay()
  {
    return edges_level_Delaunay.one_step(null_visitor);
  }

  bool try_one_step_conforming_Gabriel()
  {
    return edges_level_Gabriel.one_step(null_visitor);
  }

  bool step_by_step_conforming_Delaunay()
  {
    return edges_level_Delaunay.try_to_insert_one_point(null_visitor);
  }

  bool step_by_step_conforming_Gabriel()
  {
    return edges_level_Gabriel.try_to_insert_one_point(null_visitor);
  }
  //@}

}; // end Triangulation_conformer_2


// --- GLOBAL FUNCTIONS ---

template <class Tr>
void
make_conforming_Gabriel_2(Tr& t)
{
  typedef Triangulation_conformer_2<Tr> Conform;

  Conform conform(t);
  conform.make_conforming_Gabriel();
}

template <class Tr>
void
make_conforming_Delaunay_2(Tr& t)
{
  typedef Triangulation_conformer_2<Tr> Conform;

  Conform conform(t);
  conform.make_conforming_Delaunay();
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TRIANGULATION_CONFORMER_2_H
