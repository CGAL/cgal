// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_REFINE_EDGES_WITH_CLUSTERS_H
#define CGAL_MESH_2_REFINE_EDGES_WITH_CLUSTERS_H

#include <CGAL/Mesh_2/Refine_edges.h>
#include <CGAL/Mesh_2/Clusters.h>

namespace CGAL {

namespace Mesh_2 {
/**
 * This class is the base for the first level of Mesh_2: the edge
 * conforming level. It \e does handle clusters.
 * To handle clusters, an helping \c Clusters object is used.
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
class Refine_edges_base_with_clusters : 
    public Refine_edges_base<Tr, Is_locally_conform, Container>
{
  typedef Refine_edges_base<Tr, Is_locally_conform, Container> Super;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Point Point;
  typedef typename Tr::Geom_traits Geom_traits;

  typedef typename Geom_traits::FT FT;

  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Face_circulator Face_circulator;
  
  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

  typedef typename Clusters<Tr>::Cluster Cluster;
  typedef typename Clusters<Tr>::iterator clusters_iterator;


  /* --- private data --- */
  Clusters<Tr>& clusters;

  bool va_has_a_cluster, vb_has_a_cluster;
  bool cluster_splitted;
  Cluster ca, cb;
  clusters_iterator ca_it, cb_it;

  using Super::triangulation_ref_impl;

public:
  /** \name CONSTRUCTORS */

  Refine_edges_base_with_clusters(Tr& tr_, Clusters<Tr>& c_) 
    : Super(tr_), clusters(c_)
  {
  }


  /** \name FUNCTIONS NEEDED BY Mesher_level OVERIDDEN BY THIS CLASS. */

  Point refinement_point_impl(const Edge& edge)
  {
    typename Geom_traits::Construct_midpoint_2
      midpoint = this->triangulation_ref_impl().geom_traits().construct_midpoint_2_object();

    this->va = edge.first->vertex(Tr::cw (edge.second));
    this->vb = edge.first->vertex(Tr::ccw(edge.second));

//     std::cerr << "refinement_point_impl\n" << this->va->point() << " / "
//               << this->vb->point() << std::endl;

    va_has_a_cluster = false;
    vb_has_a_cluster = false;
    cluster_splitted = false;
    
    // true bellow to remove ca and cb because they will
    // be restored by update_cluster(...).
    if( clusters.get_cluster(this->va,this->vb,ca,ca_it) ) {
      if( clusters.get_cluster(this->vb,this->va,cb,cb_it) )
        { // both ends are clusters
          va_has_a_cluster = true;
          vb_has_a_cluster = true;
          return midpoint(this->va->point(), this->vb->point());
        }
      else {
        // va only is a cluster
        va_has_a_cluster = true;
        return split_cluster_point(this->va,this->vb,ca);
      }
    } else
    if( clusters.get_cluster(this->vb,this->va,cb,cb_it) ){
      // vb only is a cluster
      vb_has_a_cluster = true;
      return split_cluster_point(this->vb,this->va,cb);
    }else{
      // no cluster
      return midpoint(this->va->point(), this->vb->point());
    }
  };

  void after_insertion_impl(const Vertex_handle& v)
  {
#ifdef CGAL_MESH_2_DEBUG_CLUSTERS    
    std::cerr << "update_clusters" << std::endl;
    std::cerr << "va_has_a_cluster=" << va_has_a_cluster
              << std::endl
              << "vb_has_a_cluster=" << vb_has_a_cluster
              << std::endl;
    std::cerr << "clusters.size()=" << clusters.size() << std::endl;
#endif // CGAL_MESH_2_DEBUG_CLUSTERS
    Super::after_insertion_impl(v);
    if( va_has_a_cluster ) 
      clusters.update_cluster(ca,ca_it,this->va,this->vb,v,cluster_splitted);
    if( vb_has_a_cluster )
      clusters.update_cluster(cb,cb_it,this->vb,this->va,v,cluster_splitted);
#ifdef CGAL_MESH_2_DEBUG_CLUSTERS
    std::cerr << "clusters.size() after update_cluster=" 
	      << clusters.size() << std::endl;
#endif // CGAL_MESH_2_DEBUG_CLUSTERS
  }

  /**
   * Test if the edges of the boundary are locally conforming.
   * Push which that are not in the list of edges to be conformed.
   */
  Mesher_level_conflict_status
  test_point_conflict_from_superior_impl(const Point& p,
					 Zone& z)
  {
    Mesher_level_conflict_status status = NO_CONFLICT;

    // p is the circumcener of zone.parent_face.
    const FT& sq_r_of_p_parent = shortest_edge_squared_length(z.parent_face);

    for(typename Zone::Edges_iterator eit = z.boundary_edges.begin();
        eit != z.boundary_edges.end(); ++eit)
      {
        const Face_handle& fh = eit->first;
        const int& i = eit->second;

        if(fh->is_constrained(i) && 
           !this->is_locally_conform(this->tr, fh, i, p))
          {
	    const Vertex_handle& v1 = fh->vertex( this->tr.cw (i));
	    const Vertex_handle& v2 = fh->vertex( this->tr.ccw(i));

	    status = CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;

            bool v1_has_a_cluster = clusters.get_cluster(v1,v2,ca,ca_it);
            bool v2_has_a_cluster = clusters.get_cluster(v2,v1,cb,cb_it);

          if( ( v1_has_a_cluster && v2_has_a_cluster) ||
              (!v1_has_a_cluster && !v2_has_a_cluster) )
            {
              // two clusters or no cluster
              this->add_constrained_edge_to_be_conformed(v1, v2);
	      //	      status = CONF//CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
            }
          else
            {
              // only one cluster: c or c2
              if(v2_has_a_cluster)
                ca = cb;
// What Shewchuk says:
// - If the cluster is not reduced (all segments don't have the same
// length as [v1,v2]), then split the edge
// - Else, let rmin be the minimum insertion radius introduced by the
// potential split, let T be the triangle whose circumcenter
// encroaches [v1,v2] and let rg be the length of the shortest edge
// of T. If rmin >= rg, then split the edge.

              if( this->imperatively || !ca.is_reduced() ||
                  ca.rmin >=  sq_r_of_p_parent)
		this->add_constrained_edge_to_be_conformed(v1,v2);
	      else
		status = CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
            }
          }
      }  // after here edges encroached by p are in the list of edges to
         // be conformed.

    return status;
  }

private:
  /** \name Auxiliary functions to handle clusters. */

  FT shortest_edge_squared_length(Face_handle f)
  {
    typename Geom_traits::Compute_squared_distance_2 squared_distance =
      this->tr.geom_traits().compute_squared_distance_2_object();

    const Point& pa = (f->vertex(0))->point();
    const Point& pb = (f->vertex(1))->point();
    const Point& pc = (f->vertex(2))->point();

    FT a, b, c;
    a = squared_distance(pb, pc);
    b = squared_distance(pc, pa);
    c = squared_distance(pa, pb);

    return ((CGAL::min)(a, (CGAL::min)(b, c)));
  }

  Point split_cluster_point(Vertex_handle va, Vertex_handle vb, Cluster& c)
  {
    typename Geom_traits::Construct_vector_2 vector =
      this->tr.geom_traits().construct_vector_2_object();
    typename Geom_traits::Construct_scaled_vector_2 scaled_vector =
      this->tr.geom_traits().construct_scaled_vector_2_object();
    typename Geom_traits::Compute_squared_distance_2 squared_distance =
      this->tr.geom_traits().compute_squared_distance_2_object();
    typename Geom_traits::Construct_midpoint_2 midpoint =
      this->tr.geom_traits().construct_midpoint_2_object();
    typename Geom_traits::Construct_translated_point_2 translate =
      this->tr.geom_traits().construct_translated_point_2_object();

    typedef typename Geom_traits::FT FT;
    

    cluster_splitted = true;

    const Point& a = va->point();
    const Point& b = vb->point();

    if( c.is_reduced() )
      return midpoint(a, b);
    else
      {
        const Point m = midpoint(a, b);

        typename Geom_traits::Vector_2 v = vector(a,m);
        v = scaled_vector(v,CGAL_NTS sqrt(c.minimum_squared_length /
                                      squared_distance(a,b)));

        Point i = translate(a,v), i2(i);

        do {
          i = translate(a,v);
          v = scaled_vector(v,FT(2));
          i2 = translate(a,v);
        } while(squared_distance(a,i2) <= squared_distance(a,m));
        if( squared_distance(i,m) > squared_distance(m,i2) )
          i = i2;
        //here i is the best point for splitting
        return i;
      }
  }

protected:
}; // end class Refine_edges_base_with_clusters

template <
  typename Tr,
  typename Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
  typename Base = Refine_edges_base_with_clusters<Tr, 
                                                  Is_locally_conform>
>
struct Refine_edges_with_clusters : 
  public Base, 
  public details::Refine_edges_types<Tr, 
    Refine_edges_with_clusters<Tr, Is_locally_conform, Base>
    >::Edges_mesher_level
{
  typedef Refine_edges_with_clusters<Tr, Is_locally_conform, Base> Self;
  typedef typename details::Refine_edges_types<Tr,
					       Self>::Edges_mesher_level
                                 Mesher;
public:
  Refine_edges_with_clusters(Tr& t,
                             Clusters<Tr>& c,
                             Null_mesher_level& null_level)
    : Base(t, c), Mesher(null_level)
  {
  }
}; // end Refine_edges_with_clusters

} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_EDGES_WITH_CLUSTERS_H
