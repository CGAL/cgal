// Copyright (c) 2014 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot


#ifndef CGAL_STITCH_POLYGON_MESH_H
#define CGAL_STITCH_POLYGON_MESH_H

#include <CGAL/Modifier_base.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/boost/graph/helpers.h>

#include <vector>
#include <set>

namespace CGAL{

namespace Polygon_mesh_processing{

namespace internal{

template <class PolygonMesh>
struct Less_for_halfedge
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
    halfedge_descriptor;
  typedef typename boost::property_map<PolygonMesh,
                             boost::vertex_point_t>::type   Ppmap;
  typedef typename PolygonMesh::Point Point;

  Less_for_halfedge(const PolygonMesh& pmesh_,
                    const Ppmap& ppmap_)
    : pmesh(pmesh_),
      ppmap(ppmap_)
  {}

  bool operator()(halfedge_descriptor h1,
                  halfedge_descriptor h2) const
  {
    const Point& s1 = ppmap[target(opposite(h1, pmesh), pmesh)];
    const Point& t1 = ppmap[target(h1, pmesh)];
    const Point& s2 = ppmap[target(opposite(h2, pmesh), pmesh)];
    const Point& t2 = ppmap[target(h2, pmesh)];
    return
    ( s1 < t1?  std::make_pair(s1,t1) : std::make_pair(t1, s1) )
    <
    ( s2 < t2?  std::make_pair(s2,t2) : std::make_pair(t2, s2) );
  }

  const PolygonMesh& pmesh;
  const Ppmap& ppmap;
};

template <class LessHedge, class PolygonMesh, class OutputIterator>
OutputIterator
detect_duplicated_boundary_edges
(PolygonMesh& pmesh, OutputIterator out, LessHedge less_hedge)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  pmesh.normalize_border();

  typedef std::set<halfedge_descriptor, LessHedge> Border_halfedge_set;
  Border_halfedge_set border_halfedge_set(less_hedge);
  BOOST_FOREACH(halfedge_descriptor he, halfedges(pmesh))
  {
    if ( !CGAL::is_border(he, pmesh) )
      continue;
    typename Border_halfedge_set::iterator set_it;
    bool insertion_ok;
    CGAL::cpp11::tie(set_it, insertion_ok)
      = border_halfedge_set.insert(he);

    if ( !insertion_ok )
      *out++ = std::make_pair(*set_it, he);
  }
  return out;
}

template <class PolygonMesh>
struct Naive_border_stitching_modifier:
  CGAL::Modifier_base<typename PolygonMesh::HalfedgeDS>
{
  typedef typename PolygonMesh::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Vertex_handle Vertex_handle;
  typedef typename HDS::Halfedge::Base HBase;
  typedef typename PolygonMesh::Vertex::Point Point_3;

  std::vector <std::pair<Halfedge_handle,Halfedge_handle> >& hedge_pairs_to_stitch;

  Naive_border_stitching_modifier(
    std::vector< std::pair<Halfedge_handle,Halfedge_handle> >&
      hedge_pairs_to_stitch_) :hedge_pairs_to_stitch(hedge_pairs_to_stitch_)
  {}

  void update_target_vertex(Halfedge_handle h,
                            Vertex_handle v_kept,
                            CGAL::HalfedgeDS_decorator<HDS>& decorator)
  {
    Halfedge_handle start=h;
    do{
      decorator.set_vertex(h, v_kept);
      h=h->next()->opposite();
    } while( h!=start );
  }


  void operator() (HDS& hds)
  {
    std::size_t nb_hedges=hedge_pairs_to_stitch.size();

    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);

    /// Merge the vertices
    std::vector<Vertex_handle> vertices_to_delete;
    // since there might be several vertices with identical point
    // we use the following map to choose one vertex per point
    std::map<Point_3, Vertex_handle> vertices_kept;
    for (std::size_t k=0; k<nb_hedges; ++k)
    {
      Halfedge_handle h1=hedge_pairs_to_stitch[k].first;
      Halfedge_handle h2=hedge_pairs_to_stitch[k].second;

      CGAL_assertion( h1->is_border() );
      CGAL_assertion( h2->is_border() );
      CGAL_assertion( !h1->opposite()->is_border() );
      CGAL_assertion( !h2->opposite()->is_border() );

      Vertex_handle h1_tgt=h1->vertex();
      Vertex_handle h2_src=h2->opposite()->vertex();

      //update vertex pointers: target of h1 vs source of h2
      Vertex_handle v_to_keep=h1_tgt;
      std::pair<
        typename std::map<Point_3, Vertex_handle>::iterator,
        bool > insert_res =
          vertices_kept.insert( std::make_pair(v_to_keep->point(), v_to_keep) );
      if (!insert_res.second && v_to_keep!=insert_res.first->second)
      {
        v_to_keep=insert_res.first->second;
        //we remove h1->vertex()
        vertices_to_delete.push_back( h1_tgt );
        update_target_vertex(h1, v_to_keep, decorator);
      }
      if (v_to_keep!=h2_src)
      {
        //we remove h2->opposite()->vertex()
        vertices_to_delete.push_back( h2_src );
        update_target_vertex(h2->opposite(), v_to_keep, decorator);
      }
      decorator.set_vertex_halfedge(v_to_keep, h1);

      Vertex_handle h1_src=h1->opposite()->vertex();
      Vertex_handle h2_tgt=h2->vertex();

      //update vertex pointers: target of h1 vs source of h2
      v_to_keep=h2_tgt;
      insert_res =
          vertices_kept.insert( std::make_pair(v_to_keep->point(), v_to_keep) );
      if (!insert_res.second && v_to_keep!=insert_res.first->second)
      {
        v_to_keep=insert_res.first->second;
        //we remove h2->vertex()
        vertices_to_delete.push_back( h2_tgt );
        update_target_vertex(h2, v_to_keep, decorator);
      }
      if (v_to_keep!=h1_src)
      {
        //we remove h1->opposite()->vertex()
        vertices_to_delete.push_back( h1_src );
        update_target_vertex(h1->opposite(), v_to_keep, decorator);
      }
      decorator.set_vertex_halfedge(v_to_keep, h1->opposite());
    }

    /// Update next/prev of neighbor halfedges (that are not set for stiching)
    /// _______   _______
    ///        | |
    ///        | |
    /// In order to avoid having to maintain a set with halfedges to stitch
    /// we do on purpose next-prev linking that might not be useful but that
    /// is harmless and still less expensive than doing queries in a set
    for (std::size_t k=0; k<nb_hedges; ++k)
    {
      Halfedge_handle h1=hedge_pairs_to_stitch[k].first;
      Halfedge_handle h2=hedge_pairs_to_stitch[k].second;

      //link h2->prev() to h1->next()
      Halfedge_handle prev=h2->prev();
      Halfedge_handle next=h1->next();
      prev->HBase::set_next(next);
      decorator.set_prev(next, prev);

      //link h1->prev() to h2->next()
      prev=h1->prev();
      next=h2->next();
      prev->HBase::set_next(next);
      decorator.set_prev(next, prev);
    }

    /// update HDS connectivity, removing the second halfedge
    /// of each the pair and its opposite
    for (std::size_t k=0; k<nb_hedges; ++k)
    {
      Halfedge_handle h1=hedge_pairs_to_stitch[k].first;
      Halfedge_handle h2=hedge_pairs_to_stitch[k].second;

    ///Set face-halfedge relationship
      //h2 and its opposite will be removed
      decorator.set_face(h1,h2->opposite()->face());
      decorator.set_face_halfedge(h1->face(),h1);
      //update next/prev pointers
      Halfedge_handle tmp=h2->opposite()->prev();
      tmp->HBase::set_next(h1);
      decorator.set_prev(h1,tmp);
      tmp=h2->opposite()->next();
      h1->HBase::set_next(tmp);
      decorator.set_prev(tmp,h1);

    /// remove the extra halfedges
      hds.edges_erase( h2 );
    }

    //remove the extra vertices
    for(typename std::vector<Vertex_handle>::iterator
          itv=vertices_to_delete.begin(),itv_end=vertices_to_delete.end();
          itv!=itv_end; ++itv)
    {
      hds.vertices_erase(*itv);
    }
  }
};

} //end of namespace internal



/// \ingroup polyhedron_stitching_grp
/// Stitches together border halfedges in a polyhedron.
/// The halfedge to be stitched are provided in `hedge_pairs_to_stitch`.
/// Foreach pair `p` in this vector, p.second and its opposite will be removed
/// from `P`.
/// The vertices that get removed from `P` are selected as follow:
/// The pair of halfedges in `hedge_pairs_to_stitch` are processed linearly.
/// Let `p` be such a pair.
/// If the target of p.first has not been marked for deletion,
/// then the source of p.second is.
/// If the target of p.second has not been marked for deletion,
/// then the source of p.first is.
/// @tparam PolygonMesh a model of `MutableFaceGraph` and `FaceListGraph`
template <class PolygonMesh>
void stitch_borders(
  PolygonMesh& pmesh,
  std::vector <std::pair<
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor> >& 
     hedge_pairs_to_stitch)
{
  internal::Naive_border_stitching_modifier<PolygonMesh>
    modifier(hedge_pairs_to_stitch);
  pmesh.delegate(modifier);
}

/// \ingroup polyhedron_stitching_grp
/// Same as above but the pair of halfedges to be stitched are found
/// using `less_hedge`. Two halfedges `h1` and `h2` are set to be stitched
/// if `less_hedge(h1,h2)=less_hedge(h2,h1)=true`.
/// `LessHedge` is a key comparison function that is used to sort halfedges
template <class PolygonMesh, class LessHedge>
void stitch_borders(PolygonMesh& pmesh, LessHedge less_hedge)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
    halfedge_descriptor;
  std::vector <std::pair<halfedge_descriptor, halfedge_descriptor> > hedge_pairs_to_stitch;
  
  internal::detect_duplicated_boundary_edges(
    pmesh, std::back_inserter(hedge_pairs_to_stitch), less_hedge);
  stitch_borders(pmesh, hedge_pairs_to_stitch);
}

/// \ingroup polyhedron_stitching_grp
/// Same as above using the source and target points of the halfedges
/// for comparision
template <class PolygonMesh>
void stitch_borders(PolygonMesh& pmesh)
{
  typename boost::property_map<PolygonMesh, boost::vertex_point_t>::type
    ppmap = get(boost::vertex_point, pmesh);
  internal::Less_for_halfedge<PolygonMesh> less_hedge(pmesh, ppmap);
  stitch_borders(pmesh, less_hedge);
}

} //end of namespace Polygon_mesh_processing

} //end of namespace CGAL


#endif //CGAL_STITCH_POLYGON_MESH_H
