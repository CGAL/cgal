// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Ilker O. Yaz and Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_SLICER_3_H
#define CGAL_POLYGON_MESH_SLICER_3_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/tuple.h>

#include <vector>
#include <set>

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/variant.hpp>

#include <CGAL/boost/graph/split_graph_into_polylines.h>

namespace CGAL {

/// \ingroup PkgPolygonMeshProcessing
/// Function object that can compute the intersection of planes with
/// a triangulated surface mesh.
/// \tparam TriangleMesh must be a model of `FaceGraph` and `HalfedgeListGraph`
/// \tparam Traits must be a model of `AABBGeomTraits`
/// \tparam VertexPointPmap is a model of `ReadablePropertyMap` with
///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
///         `Traits::Point_3` as value type
/// \tparam AABBTree must be an instanciation of `CGAL::AABB_tree` able to handle
///         the edges of TriangleMesh, having its `edge_descriptor` as primitive id.
/// Depends on \ref PkgAABB_treeSummary
/// \todo use a custom traversal traits
/// \todo use dedicated predicates if the plane is axis aligned
template<class TriangleMesh,
  class Traits,
  class VertexPointPmap = typename boost::property_map< TriangleMesh, vertex_point_t>::type,
  class AABBTree = AABB_tree<
                       AABB_traits<Traits,
                         AABB_halfedge_graph_segment_primitive<TriangleMesh> > > >
class Polygon_mesh_slicer_3
{
/// Polygon_mesh typedefs
  typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;

/// Geometric typedefs
  typedef typename Traits::Plane_3                                      Plane_3;
  typedef typename Traits::Segment_3                                  Segment_3;
  typedef typename Traits::Intersect_3                              Intersect_3;
  typedef typename Traits::Point_3                                      Point_3;

  enum Intersection_type_enum { POINT, SRC_VERTEX, SEGMENT };
  typedef std::pair<edge_descriptor, Intersection_type_enum>  Intersection_type;

/// typedefs for internal graph to get connectivity of the polylines
  typedef boost::variant<vertex_descriptor, edge_descriptor>     AL_vertex_info;
  typedef boost::adjacency_list <
                              boost::vecS,
                              boost::vecS,
                              boost::undirectedS,
                              AL_vertex_info >                         AL_graph;
  typedef typename AL_graph::vertex_descriptor             AL_vertex_descriptor;
  typedef std::pair<AL_vertex_descriptor, AL_vertex_descriptor>  AL_vertex_pair;

  // compare the faces using the halfedge descriptors
  struct Compare_face{
    TriangleMesh& m_tmesh;
    Compare_face(TriangleMesh& tmesh)
      :m_tmesh(tmesh)
    {}

    bool operator()(halfedge_descriptor hd1, halfedge_descriptor hd2) const
    {
      return face(hd1,m_tmesh) < face(hd2,m_tmesh);
    }
  };

  typedef std::map< halfedge_descriptor, AL_vertex_pair, Compare_face > AL_edge_map;

  template <class OutputIterator>
  struct Polyline_visitor{
    AL_graph& al_graph;
    TriangleMesh& m_tmesh;
    const Plane_3& m_plane;
    VertexPointPmap m_vpmap;
    const Traits& m_traits;
    OutputIterator out;

    Polyline_visitor( TriangleMesh& tmesh,
                      AL_graph& al_graph,
                      const Plane_3& plane,
                      VertexPointPmap vpmap,
                      const Traits& traits,
                      OutputIterator out)
      : al_graph(al_graph)
      , m_tmesh(tmesh)
      , m_plane(plane)
      , m_vpmap(vpmap)
      , m_traits(traits)
      , out(out)
    {}

    std::vector< Point_3 > current_poly;
    void start_new_polyline()
    {
      current_poly.clear();
    }
    void add_node(AL_vertex_descriptor node_id)
    {
      AL_vertex_info v = al_graph[node_id];
      if (const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&v) )
      {
        current_poly.push_back( get(m_vpmap, *vd_ptr) );
      }
      else
      {
        edge_descriptor ed = boost::get<edge_descriptor>(v);
        Segment_3 s(
          get(m_vpmap, source(ed, m_tmesh)),
          get(m_vpmap,target(ed, m_tmesh))
        );
        Intersect_3 intersection = m_traits.intersect_3_object();
        typename cpp11::result_of<Intersect_3(Plane_3, Segment_3)>::type
          inter = intersection(m_plane, s);
        CGAL_assertion( inter );
        const Point_3* pt_ptr = boost::get<Point_3>(&(*inter));
        current_poly.push_back( *pt_ptr );
      }
    }

    void end_polyline()
    {
      CGAL_assertion(!current_poly.empty());
      *out++=current_poly;
    }
  };

/// member variables
  const AABBTree* m_tree_ptr;
  TriangleMesh& m_tmesh;
  VertexPointPmap m_vpmap;
  Traits m_traits;
  bool m_own_tree;

/// Convenience graph functions
  edge_descriptor opposite_edge(edge_descriptor ed) const
  {
    return edge( opposite( halfedge(ed, m_tmesh), m_tmesh), m_tmesh );
  }

  edge_descriptor next_edge(edge_descriptor ed) const
  {
    return edge( next( halfedge(ed, m_tmesh), m_tmesh), m_tmesh );
  }

  edge_descriptor next_of_opposite_edge(edge_descriptor ed) const
  {
    return edge( next( opposite( halfedge(ed, m_tmesh), m_tmesh), m_tmesh), m_tmesh );
  }

  face_descriptor opposite_face(edge_descriptor ed) const
  {
    return face( opposite( halfedge(ed, m_tmesh), m_tmesh), m_tmesh);
  }

  /// given an edge intersected by `plane` indicates whether the edge is
  /// intersected in its interior, at an endpoint or if it's coplanar
  Intersection_type
  classify_edge(edge_descriptor ed, const Plane_3& plane) const
  {
    typename Traits::Oriented_side_3 oriented_side = m_traits.oriented_side_3_object();
    Oriented_side src = oriented_side(plane, get(m_vpmap, source(ed,m_tmesh)) );
    Oriented_side tgt = oriented_side(plane, get(m_vpmap, target(ed,m_tmesh)) );

    if (src==ON_ORIENTED_BOUNDARY)
    {
      return tgt==ON_ORIENTED_BOUNDARY?
                  Intersection_type(ed, SEGMENT):
                  Intersection_type(ed, SRC_VERTEX);
    }

    if (tgt==ON_ORIENTED_BOUNDARY)
      return Intersection_type(opposite_edge(ed), SRC_VERTEX);

    CGAL_assertion (src!=tgt);
    return Intersection_type(ed,POINT);
  }

  /// handle edge insertion in the adjacency_list graph
  /// we add an edge betweem two edge_descriptor if they
  /// share a common facet
  void update_al_graph_connectivity(
    edge_descriptor ed,
    AL_vertex_descriptor vd,
    AL_edge_map& al_edge_map,
    AL_graph& al_graph) const
  {
      typename AL_edge_map::iterator itm;
      bool new_insertion;

      halfedge_descriptor hd=halfedge(ed, m_tmesh);

      if (face(hd, m_tmesh)!=graph_traits::null_face())
      {
        cpp11::tie(itm, new_insertion) =
          al_edge_map.insert( std::pair< halfedge_descriptor, AL_vertex_pair >
            (hd, AL_vertex_pair(vd, AL_graph::null_vertex())) );
        if (!new_insertion)
        {
          CGAL_assertion(itm->second.second==AL_graph::null_vertex());
          itm->second.second=vd;
          add_edge( itm->second.first,
                    itm->second.second,
                    al_graph);
        }
      }

      hd=opposite(hd, m_tmesh);
      if (face(hd, m_tmesh)!=graph_traits::null_face())
      {
        cpp11::tie(itm, new_insertion) =
          al_edge_map.insert( std::pair< halfedge_descriptor, AL_vertex_pair >
            (hd, AL_vertex_pair(vd, AL_graph::null_vertex())) );
        if (!new_insertion)
        {
          CGAL_assertion(itm->second.second==AL_graph::null_vertex());
          itm->second.second=vd;
          add_edge( itm->second.first,
                    itm->second.second,
                    al_graph);
        }
      }
  }

public:

  /**
  * Constructor using all `edges(tmesh)` to initialize the
  * internal `AABB_tree`.
  * @param tmesh the triangulated surface mesh to be cut.
  *              it must be valid and non modified as long
  *              as the functor is used
  * @param vpmap an intance of the vertex point property map
  * @param traits a traits class instance
  */
  Polygon_mesh_slicer_3(const TriangleMesh& tmesh,
                        VertexPointPmap vpmap,
                        const Traits& traits = Traits())
  : m_tmesh(const_cast<TriangleMesh&>(tmesh))
  , m_vpmap(vpmap)
  , m_traits(traits)
  , m_own_tree(true)
  {
    m_tree_ptr = new AABBTree(edges(m_tmesh).first,
                              edges(m_tmesh).second,
                              m_tmesh,
                              m_vpmap);
  }

  /**
  * Constructor using an `AABB_tree` provided by the user.
  * @param tmesh the triangulated surface mesh to be cut.
  *              it must be valid and non modified as long
  *              as the functor is used
  * @param tree must be initialized with all the edge of `tmesh`
  * @param vpmap an intance of the vertex point property map
  * @param traits a traits class instance
  */
  Polygon_mesh_slicer_3(const TriangleMesh& tmesh,
                        const AABBTree& tree,
                        VertexPointPmap vpmap,
                        const Traits& traits = Traits())
    : m_tree_ptr(&tree)
    , m_tmesh(const_cast<TriangleMesh&>(tmesh))
    , m_vpmap(vpmap)
    , m_traits(traits)
    , m_own_tree(false)
  { }

  /**
  * Constructor using all `edges(tmesh)` to initialize the
  * internal `AABB_tree`. The vertex point property map used
  * is `get(boost::vertex_point, tmesh)`
  * @param tmesh the triangulated surface mesh to be cut.
  *              it must be valid and non modified as long
  *              as the functor is used
  * @param traits a traits class instance
  */
  Polygon_mesh_slicer_3(const TriangleMesh& tmesh,
                        const Traits& traits = Traits())
  : m_tmesh(const_cast<TriangleMesh&>(tmesh))
  , m_vpmap(get(boost::vertex_point, m_tmesh))
  , m_traits(traits)
  , m_own_tree(true)
  {
    m_tree_ptr = new AABBTree(edges(m_tmesh).first,
                              edges(m_tmesh).second,
                              m_tmesh,
                              m_vpmap);
  }

  /**
  * Constructor using an `AABB_tree` provided by the user.
  * The vertex point property map used is `get(boost::vertex_point, tmesh)`
  * @param tmesh the triangulated surface mesh to be cut.
  *              it must be valid and non modified as long
  *              as the functor is used
  * @param tree must be initialized with all the edge of `tmesh`
  * @param traits a traits class instance
  */
  Polygon_mesh_slicer_3(const TriangleMesh& tmesh,
                        const AABBTree& tree,
                        const Traits& traits = Traits())
    : m_tree_ptr(&tree)
    , m_tmesh(const_cast<TriangleMesh&>(tmesh))
    , m_vpmap(get(boost::vertex_point, m_tmesh))
    , m_traits(traits)
    , m_own_tree(false)
  { }

  /**
   * Construct the intersecting polylines of `plane` with the input triangulated surface mesh.
   * @tparam OutputIterator an output iterator accepting polylines. A polyline is provided as `std::vector<Kernel::Point_3>`.
   *                        A polyline is closed if its first and last points are identical.
   * @param plane the plane to intersect the triangulated surface mesh with
   * @param out output iterator of polylines
   */
  template <class OutputIterator>
  OutputIterator operator() (const Plane_3& plane,
                             OutputIterator out) const
  {
    CGAL_precondition(!plane.is_degenerate());

    /// get all edges intersected by the plane
    std::vector<edge_descriptor> intersected_edges;
    m_tree_ptr->all_intersected_primitives(plane, std::back_inserter(intersected_edges));

    std::set<edge_descriptor> all_coplanar_edges;
    std::vector<edge_descriptor> iedges;
    typedef std::map<vertex_descriptor, AL_vertex_descriptor> Src_vertices_map;
    Src_vertices_map src_vertices;

    AL_graph al_graph; //< output graph

    // classify each intersected edge according to how the plane intersects it
    BOOST_FOREACH(edge_descriptor ed, intersected_edges)
    {
      Intersection_type itype=classify_edge(ed, plane);
      switch(itype.second)
      {
        case SEGMENT:
          all_coplanar_edges.insert(ed);
        break;
        case SRC_VERTEX:
        {
          //insert vertices and make sure it is a new one
          typename Src_vertices_map::iterator it_insert;
          bool is_new;
          cpp11::tie(it_insert, is_new) =
            src_vertices.insert(
              std::pair<vertex_descriptor,AL_vertex_descriptor>(
                source(itype.first,m_tmesh), AL_graph::null_vertex()
              )
            );
          if (is_new)
          {
            it_insert->second=add_vertex(al_graph);
            al_graph[it_insert->second]=it_insert->first;
          }
        }
        break;
        default:
          iedges.push_back(ed);
      }
    }

    Compare_face less_face(m_tmesh);
    AL_edge_map al_edge_map( less_face );

    /// Filter coplanar edges: we consider only coplanar edges incident to one non-coplanar facet
    ///   for each such edge, add the corresponding nodes in the adjacency-list graph as well as
    ///   the edge
    BOOST_FOREACH(const edge_descriptor ed, all_coplanar_edges)
    {
      if (  face(halfedge(ed, m_tmesh), m_tmesh)==graph_traits::null_face() ||
            opposite_face(ed)==graph_traits::null_face()  ||
            !all_coplanar_edges.count( next_edge(ed) ) ||
            !all_coplanar_edges.count( next_of_opposite_edge(ed) ) )
      {
        typename Src_vertices_map::iterator it_insert1, it_insert2;
        bool is_new;

        /// Each coplanar edge is connecting two nodes
        //  handle source
        cpp11::tie(it_insert1, is_new) =
          src_vertices.insert(
              std::pair<vertex_descriptor,AL_vertex_descriptor>(
                source(ed,m_tmesh), AL_graph::null_vertex()
              )
          );
        if (is_new)
        {
          it_insert1->second=add_vertex(al_graph);
          al_graph[it_insert1->second]=it_insert1->first;
        }
        //  handle target
        cpp11::tie(it_insert2, is_new) =
          src_vertices.insert(
              std::pair<vertex_descriptor,AL_vertex_descriptor>(
                target(ed,m_tmesh), AL_graph::null_vertex()
              )
          );
        if (is_new)
        {
          it_insert2->second=add_vertex(al_graph);
          al_graph[it_insert2->second]=it_insert2->first;
        }
        /// add the edge into the adjacency-list graph
        CGAL_assertion( it_insert1->second!=AL_graph::null_vertex() );
        CGAL_assertion( it_insert2->second!=AL_graph::null_vertex() );
        add_edge(it_insert1->second, it_insert2->second, al_graph);
      }
    }

    /// for each edge intersected in its interior, creates a node in
    /// an adjacency-list graph and put an edge between two such nodes
    /// when the corresponding edges shares a common face
    BOOST_FOREACH(edge_descriptor ed, iedges)
    {
      AL_vertex_descriptor vd=add_vertex(al_graph);
      al_graph[vd]=ed;
      update_al_graph_connectivity(ed, vd, al_edge_map, al_graph);
    }

    /// If one of the node above is not connected in its two incident faces
    /// then it must be connected to a vertex (including those in the set
    /// of coplanar edges)
    typedef std::pair<halfedge_descriptor, AL_vertex_pair> Halfedge_and_vertices;
    BOOST_FOREACH(Halfedge_and_vertices hnv,al_edge_map)
    {
      if (hnv.second.second==AL_graph::null_vertex())
      {
        //get the edge and test opposite vertices (if the edge is not on the boundary)
        vertex_descriptor vd = target( next(hnv.first, m_tmesh), m_tmesh);
        typename Src_vertices_map::iterator itv=src_vertices.find(vd);
        CGAL_assertion( itv!=src_vertices.end() );
        add_edge(itv->second, hnv.second.first, al_graph);
      }
    }

    CGAL_assertion(num_vertices(al_graph)==iedges.size()+src_vertices.size());

    /// now assemble the edges of al_graph to define polylines,
    /// putting them in the output iterator
    Polyline_visitor<OutputIterator> visitor(m_tmesh, al_graph, plane, m_vpmap, m_traits, out);
    split_graph_into_polylines(al_graph, visitor);
    return visitor.out;
  }

  ~Polygon_mesh_slicer_3()
  {
    if (m_own_tree) delete m_tree_ptr;
  }
};

}// end of namespace CGAL
#endif //CGAL_POLYGON_MESH_SLICER_3_H
