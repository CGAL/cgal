// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// Copyright (c) 2010       GeometryFactory Sarl (France)
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_2_CLUSTERS_H
#define CGAL_MESH_2_CLUSTERS_H

#include <CGAL/Filter_circulator.h>
#include <CGAL/Unique_hash_map.h>

#include <utility>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

namespace Mesh_2 
{

  namespace details 
  {
    template <class Tr>
    class Is_edge_constrained {
      const Tr* tr_;
    public:
      typedef Is_edge_constrained<Tr> Self;
      typedef typename Tr::Edge_circulator Edge_circulator;
      
      Is_edge_constrained(const Tr& tr) : tr_(&tr)
      {}

      bool operator()(const Edge_circulator& ec) const
      {
        return tr_->is_constrained(*ec);
      }
    };
  } // end namespace details

template <class Tr>
class Clusters
{
  typedef typename Tr::Vertex_handle          Vertex_handle;
  typedef typename Tr::Point                  Point;
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef FT      Squared_length; /**<This typedef is used to remind that
                                     the length is squared. */
  typedef typename Tr::Edge_circulator Edge_circulator;
  
  /**
   *  Special type: filtered circulator that returns only constrained
   *  edges.
   */
  typedef Filter_circulator<Edge_circulator,
                            details::Is_edge_constrained<Tr> >
    Constrained_edge_circulator;

public:
  /** \name Clusters public types */

  /**
   * \c Cluster register several informations about clusters.
   * A cluster is a set of vertices v_i incident to one vertice
   * v_0, so that angles between segments [v_0, v_i] is less than 60
   * degres.
   */
  struct Cluster {
    bool reduced ; /**< Is the cluster reduced? */

    /** 
     * Smallest_angle gives the two vertices defining the
     * smallest angle in the cluster.
     */
    std::pair<Vertex_handle, Vertex_handle> smallest_angle;

    FT rmin; // @fixme: rmin has no meaning if reduced=false!!!
    Squared_length minimum_squared_length;

    /**
     * The following map tells what vertices are in the cluster and if
     * the corresponding segment has been splitted once.
     */
    typedef std::map<Vertex_handle, bool> Vertices_map;
    Vertices_map vertices;

    bool is_reduced() const {
      return reduced;
    }

    bool is_reduced(const Vertex_handle v) {
      return vertices[v];
    }
  };
private:
  /** \name Clusters associated types */

  typedef std::multimap<Vertex_handle, Cluster> Cluster_map;
  typedef typename Cluster_map::value_type Cluster_map_value_type;

  template <class Pair>
  struct Pair_get_first: public std::unary_function<Pair,
                                                    typename Pair::first_type>
  {
    typedef typename Pair::first_type result;
    const result& operator()(const Pair& p) const
    {
      return p.first;
    }
  };

  typedef typename Cluster::Vertices_map Cluster_vertices_map;

private:
  /* --- protected datas --- */

  Tr& tr; /**< The triangulation itself. */

  /**
   * Multimap \c Vertex_handle -> \c Cluster
   * Each vertex can have several clusters. 
   */
  Cluster_map cluster_map;

public:
  typedef typename Cluster_map::const_iterator const_iterator;
  typedef typename Cluster_map::iterator iterator;

  Clusters(Tr& tr_) : tr(tr_)
  {
  }

  /** For all vertices, calls create_clusters_of_vertex(). */
  void create_clusters() {
    create_clusters(typename Tr::Constraint_hierarchy_tag());
  }

  // function that depends of Tr::Constraint_hierarchy_tag
  template <typename Constraint_hierarchy_tag>
  void create_clusters(Constraint_hierarchy_tag) {
    cluster_map.clear();
    for(typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end();
        vit++)
    {
      create_clusters_of_vertex(vit);
    }
  }

  void create_clusters(Tag_true) {
    cluster_map.clear();
    Unique_hash_map<Vertex_handle,bool> created(false);
    for(typename Tr::Subconstraint_iterator it = tr.subconstraints_begin();
        it != tr.subconstraints_end(); ++it) {
      Vertex_handle vh = it->first.first;
      if(!created[vh]){
        created[vh] = true;
        create_clusters_of_vertex(vh);
      }

      vh = it->first.second;
      if(!created[vh]){
        created[vh] = true;
        create_clusters_of_vertex(vh);
      }
    }
  }

private:
  /**
   * Computes clusters of the vertex \c v, using the auxiliary function
   * construct_cluster().
   */
  void create_clusters_of_vertex(const Vertex_handle v);

  /**
   * Adds the sequence [\c begin, \c end] to the cluster \c c and adds it 
   * to the clusters of the vertex \c v.
   */
  void construct_cluster(const Vertex_handle v,
                         const Constrained_edge_circulator& begin,
                         const Constrained_edge_circulator& end,
                         Cluster c = Cluster());

public:
  /** \name Functions to manage clusters during the refinement process. */

  /** 
   * Update the cluster of [\c va,\c vb], putting \c vm instead of \c vb.
   * If reduction=false, the edge [va,vm] is not set reduced. 
   */
  void update_cluster(Cluster& c, iterator it,
                      const Vertex_handle va, const Vertex_handle vb,
                      const Vertex_handle vm,
                      bool reduction = true);

  /**
   * Returns the cluster of [\c va,\c vb] in \c c and return true
   * if it is in a cluster. Returns also a const_iterator in \c it.
   */
  bool get_cluster(const Vertex_handle va, const Vertex_handle vb,
                   Cluster& c, iterator& it);

  /** Const version of get_cluster(). */
  bool get_cluster(const Vertex_handle va, const Vertex_handle vb,
                   Cluster& c, const_iterator& it) const;

  /** \name Auxiliary functions that return a boolean. */

  /**
   * Tells if the angle <pleft, pmiddle, pright> is less than 60 degres.
   * Uses squared_cosine_of_angle_times_4() and used by
   * create_clusters_of_vertex().
   */
  bool is_small_angle(const Point& pleft,
                      const Point& pmiddle,
                      const Point& pright) const;

private:
  /** \name Helping computing functions */

  /** Returns the squared cosine of the angle <pleft, pmiddle, pright>
      times 4. */
  FT squared_cosine_of_angle_times_4(const Point& pleft,
                                     const Point& pmiddle,
                                     const Point& pright) const;

  /** Helper functions to access the two vertices of an Edge
      source is the vertex around which the circulator turns. */
  //@{
  Vertex_handle source(const Edge_circulator& ec) const
  {
    return ec->first->vertex(tr.cw(ec->second));
  }

  Vertex_handle target(const Edge_circulator& ec) const
  {
    return ec->first->vertex(tr.ccw(ec->second));
  }
  //@}

public:
  /** \name CONST ACCESS FUNCTIONS */
  typedef typename boost::transform_iterator<
    Pair_get_first<typename Cluster_map::value_type>,
    typename Cluster_map::const_iterator>
  Cluster_vertices_iterator;

  typedef typename boost::transform_iterator<
    Pair_get_first<typename Cluster_vertices_map::value_type>,
    typename Cluster_vertices_map::const_iterator>
  Vertices_in_cluster_iterator;

  int size() const
  {
    return cluster_map.size();
  }

  Cluster_vertices_iterator clusters_vertices_begin() const
  {
    return Cluster_vertices_iterator(cluster_map.begin());
  }

  Cluster_vertices_iterator clusters_vertices_end() const
  {
    return Cluster_vertices_iterator(cluster_map.end());
  }

  unsigned int number_of_clusters_at_vertex(const Vertex_handle& vh) const 
  {
    typedef typename Cluster_map::const_iterator Iterator;
    typedef std::pair<Iterator, Iterator> Range;
    Range range = cluster_map.equal_range(vh);
    return std::distance(range.first, range.second);
  }

  // returns the sequence of vertices bellonging to the n-th cluster of vh
  std::pair<Vertices_in_cluster_iterator, Vertices_in_cluster_iterator>
  vertices_in_cluster_sequence(const Vertex_handle& vh,
                               const unsigned int n) const
  {
    typedef typename Cluster_map::const_iterator Iterator;
    typedef std::pair<Iterator, Iterator> Range;

    Range range = cluster_map.equal_range(vh);
    Iterator first = range.first;
    std::advance(first, n);
    const Cluster& c = first->second;

    return
      std::make_pair(Vertices_in_cluster_iterator(c.vertices.begin()),
                     Vertices_in_cluster_iterator(c.vertices.end()));
  }

}; // end class Clusters

template <typename Tr>
void Clusters<Tr>::
update_cluster(Cluster& c, iterator it, Vertex_handle va,
               Vertex_handle vb, Vertex_handle vm, bool reduction)
{
  typename Geom_traits::Compute_squared_distance_2 squared_distance =
    tr.geom_traits().compute_squared_distance_2_object();

  cluster_map.erase(it);

  c.vertices.erase(vb);
  c.vertices[vm] = reduction;

  if(vb==c.smallest_angle.first)
    c.smallest_angle.first = vm;
  if(vb==c.smallest_angle.second)
    c.smallest_angle.second = vm;

  FT l = squared_distance(va->point(),vm->point());
  if(l<c.minimum_squared_length)
    c.minimum_squared_length = l;

  if(!c.is_reduced())
    {
      typename Cluster::Vertices_map::iterator it = c.vertices.begin();
      while(it!=c.vertices.end() && c.is_reduced(it->first))
        ++it; // @todo: use std::find and an object class
      if(it==c.vertices.end())
        c.reduced = true;
    }

  if(c.is_reduced())
    c.rmin = squared_distance(c.smallest_angle.first->point(),
                              c.smallest_angle.second->point())/FT(4);
  cluster_map.insert(Cluster_map_value_type(va,c));
}

template <typename Tr>
bool Clusters<Tr>::
get_cluster(Vertex_handle va, Vertex_handle vb, Cluster& c,
            const_iterator& it) const
{
  typedef std::pair<const_iterator, const_iterator> Range;

  Range range = cluster_map.equal_range(va);

  for(it = range.first; it != range.second; it++)
    {
      const Cluster &cl = it->second;
      if(cl.vertices.find(vb)!=cl.vertices.end()) {
        c = it->second;
        return true;
      }
    }
  return false;
}

template <typename Tr>
bool Clusters<Tr>::
get_cluster(Vertex_handle va, Vertex_handle vb, Cluster& c,
            iterator& it) 
{
  typedef std::pair<iterator, iterator> Range;

  Range range = cluster_map.equal_range(va);

  for(it = range.first; it != range.second; it++)
    {
      const Cluster &cl = it->second;
      if(cl.vertices.find(vb)!=cl.vertices.end()) {
        c = it->second;
        return true;
      }
    }
  return false;
}


template <typename Tr>
void Clusters<Tr>::
create_clusters_of_vertex(const Vertex_handle v)
{
  details::Is_edge_constrained<Tr> test(tr);

  Constrained_edge_circulator begin(tr.incident_edges(v),test);

  // This circulator represents all constrained edges around the
  // vertex v. An edge [v,v'] is represented by the vertex v'.

  if(begin == 0) return; // if there is only one vertex

  Constrained_edge_circulator
    current(begin), next(begin), cluster_begin(begin);
  ++next; // next is always just after current.
  if(current == next) return;

  bool in_a_cluster = false;
  do
    {
      if(is_small_angle(target(current)->point(), v->point(),
                        target(next)->point()))
        {
          if(!in_a_cluster)
            {
              // at this point, current is the beginning of a cluster
              in_a_cluster = true;
              cluster_begin = current;
            }
        }
      else
        if(in_a_cluster)
          {
            // at this point, current is the end of a cluster and
            // cluster_begin is its beginning
            construct_cluster(v, cluster_begin, current);
            in_a_cluster = false;
          }
      ++next;
      ++current;
    } while( current!=begin );
  if(in_a_cluster)
    {
      Cluster c;
      iterator it;
      if(get_cluster(v, target(begin), c, it))
        {
          // get the cluster and erase it from the clusters map
          cluster_map.erase(it);
          construct_cluster(v, cluster_begin, begin, c);
        }
      else
        construct_cluster(v, cluster_begin, current);
    }
}

template <typename Tr>
void Clusters<Tr>::
construct_cluster(Vertex_handle v,
                  const Constrained_edge_circulator& begin,
                  const Constrained_edge_circulator& end,
                  Cluster c)
{
  typename Geom_traits::Compute_squared_distance_2 squared_distance =
    tr.geom_traits().compute_squared_distance_2_object();

  if(c.vertices.empty())
    {
      c.reduced = false;
      // c.rmin is not initialized because
      // reduced=false!
      c.minimum_squared_length =
        squared_distance(v->point(), target(begin)->point());
      Constrained_edge_circulator second(begin);
      ++second;
      c.smallest_angle.first = target(begin);
      c.smallest_angle.second = target(second);
    }

  const bool all_edges_in_cluster = (begin == end); // tell if all incident edges
                                              // are in the cluster
  const Point& vp = v->point();

  FT greatest_cosine =
    squared_cosine_of_angle_times_4(c.smallest_angle.first->point(),
                                    v->point(),
                                    c.smallest_angle.second->point());

  bool one_full_loop_is_needed = all_edges_in_cluster;

  bool stop = false;
  Constrained_edge_circulator circ(begin);
  Constrained_edge_circulator next(begin);
  while(!stop)
  {
    c.vertices[target(circ)] = false;
    Squared_length l = squared_distance(vp,
                                        target(circ)->point());
    c.minimum_squared_length =
      (std::min)(l,c.minimum_squared_length);

    if(circ!=end || one_full_loop_is_needed)
    {
      FT cosine =
        squared_cosine_of_angle_times_4(target(circ)->point(),
                                        v->point(),
                                        target(next)->point());
      if(cosine>greatest_cosine)
      {
        greatest_cosine = cosine;
        c.smallest_angle.first = target(circ);
        c.smallest_angle.second = target(next);
      }
    }

    if(one_full_loop_is_needed) {
      one_full_loop_is_needed = false;
    } else {
      stop = (circ == end);
    }
    ++circ;
    ++next;
  }

  typedef typename Cluster_map::value_type Value_key_pair;
  cluster_map.insert(Value_key_pair(v,c));
}

template <typename Tr>
bool Clusters<Tr>::
is_small_angle(const Point& pleft,
               const Point& pmiddle,
               const Point& pright) const
{
  typename Geom_traits::Angle_2 angle = 
    tr.geom_traits().angle_2_object();
  typename Geom_traits::Orientation_2 orient =
    tr.geom_traits().orientation_2_object();

  if( angle(pleft, pmiddle, pright)==OBTUSE )
    return false;
  if( orient(pmiddle,pleft,pright)==RIGHT_TURN)
    return false;

  FT cos_alpha = squared_cosine_of_angle_times_4(pleft, pmiddle,
                                                 pright);

  if(cos_alpha > 1)
    {
      return true; //the same cluster
    }
  else
    {
      return false; //another cluster
    }
}

template <typename Tr>
typename Clusters<Tr>::FT
Clusters<Tr>::
squared_cosine_of_angle_times_4(const Point& pb, const Point& pa,
                                const Point& pc) const
{
  typename Geom_traits::Compute_squared_distance_2 squared_distance =
    tr.geom_traits().compute_squared_distance_2_object();

  const FT
    a = squared_distance(pb, pc),
    b = squared_distance(pa, pb),
    c = squared_distance(pa, pc);

  const FT num = a-(b+c);

  return (num*num)/(b*c);
}
  
} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_CLUSTERS_H
