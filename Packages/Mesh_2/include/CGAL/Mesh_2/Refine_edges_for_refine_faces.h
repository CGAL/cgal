// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_REFINE_EDGES_FOR_REFINE_FACES_H
#define CGAL_MESH_2_REFINE_EDGES_FOR_REFINE_FACES_H

#include <CGAL/Mesh_2/Refine_edges_with_clusters.h>
#include <CGAL/Mesh_2/Clusters.h>

#include <set>
#include <list>
#include <iterator>

namespace CGAL {

namespace Mesh_2 {
/**
 * This class is the base for the first level of Mesh_2: the edge
 * conforming level.
 * This version is a modified Refine_edges_with_clusters, that adds an hack
 * used by Refine_edges_visitor.
 *
 * \param Tr is the type of triangulation on which the level acts.
 * \param Is_locally_conform defines the locally conform criterion: Gabriel
 *        or Delaunay. It defaults to the Garbriel criterion.
 * \param Container is the type of container. It defaults to a filtered
 *        queue of \c Vertex_handle pair (see \c Filtered_queue_container).
 */
template <
  class Tr,
  class Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
  class Container = 
    typename details::Refine_edges_base_types<Tr>::Default_container
>
class Refine_edges_base_for_refine_faces : 
    public Refine_edges_base_with_clusters<Tr,
                                           Is_locally_conform,
                                           Container>
{
  typedef Refine_edges_base_with_clusters<Tr,
                                          Is_locally_conform,
                                          Container> Super;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Point Point;
  typedef typename Tr::Geom_traits Geom_traits;

  typedef typename Geom_traits::FT FT;

  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Face_circulator Face_circulator;
  
  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

public:
  typedef std::set<Face_handle> Special_zone;

private: // --- private data member
  Special_zone& special_zone;

public:
  /** \name CONSTRUCTORS */

  Refine_edges_base_for_refine_faces(Tr& tr_, Clusters<Tr>& c_,
                                     Special_zone& sp_)
    : Super(tr_, c_), special_zone(sp_)
  {
  }


  /* \name FUNCTIONS NEEDED BY \c Mesher_level OVERIDDEN BY THIS CLASS. */

  /** Unmark as constrained.
      Insert flipped faces in special_zone.
      This function OVERIDES  Super::do_before_conflicts().
  */
  void do_before_conflicts(const Edge& e, const Point&)
  {
    std::cerr << e.first->vertex(Tr::cw(e.second))->point() << std::endl;
    this->tr.remove_constrained_edge(e.first, e.second,
                                     std::inserter(special_zone,
                                                   special_zone.begin()));
    std::cerr << this->va->point() << " / "
              << this->vb->point() << std::endl;
    std::cerr << e.first->vertex(Tr::cw(e.second))->point() << std::endl;
  }

  /** 
   * Finds which faces are intersected by the two sub-segments. These faces
   * are inserted in special_zone.
   * Calls Super::do_after_insertion(), that will insert the two
   * subsegments.
   */
  void do_after_insertion(const Vertex_handle& v)
  {
    typename Tr::List_edges dummy_edges;
    typedef typename Tr::List_faces List_faces;
    List_faces list_faces;
    Vertex_handle dummy_vh;

    std::cerr << this->va->point() << " / "
              << this->vb->point() << std::endl;

    CGAL_assertion_code( bool should_be_false = )
    this->tr.find_intersected_faces(this->va, this->vb,
                                    list_faces,
                                    dummy_edges,
                                    dummy_edges,
                                    dummy_vh);
    CGAL_assertion( should_be_false == false );
    for(typename List_faces::const_iterator fit = list_faces.begin();
        fit != list_faces.end();
        ++fit)
      special_zone.insert(*fit);

    Super::do_after_insertion(v);
  }

  /** Clean special_zone. */
  void do_after_no_insertion(const Edge& e, const Point& p,
                             const Zone& z)
  {
    special_zone.clear();

    Super::do_after_no_insertion(e, p, z);
  }

}; // end Refine_edges_for_refine_faces

template <
  typename Tr,
  typename Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
  typename Base = Refine_edges_base_for_refine_faces<Tr, 
                                                     Is_locally_conform>
>
struct Refine_edges_for_refine_faces : 
  public Base, 
  public details::Refine_edges_types<Tr, 
    Refine_edges_for_refine_faces<Tr, Is_locally_conform, Base> >::
    Edges_mesher_level
{
  typedef Refine_edges_for_refine_faces<Tr, Is_locally_conform, Base> Self;
  typedef typename details::Refine_edges_types<Tr, Self>::Edges_mesher_level
                              Mesher;
public:
  typedef typename Self::Special_zone Special_zone;

  Refine_edges_for_refine_faces(Tr& t,
                                Clusters<Tr>& c,
                                Special_zone& special_zone,
                                Null_mesher_level& null_level)
    : Base(t, c, special_zone), Mesher(null_level)
  {
  }
}; // end Refine_edges_for_refine_faces

} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_EDGES_FOR_REFINE_FACES_H
