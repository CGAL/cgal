// Copyright (c) 2011 GeometryFactory (France).
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
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H
#define CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/internal/corefinement/utils.h>

#include <CGAL/property_map.h>
#include <boost/optional.hpp>
#include <boost/next_prior.hpp>
#include <boost/foreach.hpp>

#include <fstream>
#include <sstream>
#include <iterator>

//  TODOCUMENT
//  --We suppose that the two input polyhedra are triangulated orientable surfaces.
//  --Any polyhedron defines two bounding volumes: one inside and one outside.
//    The convention used is the following: the normal of a triangle always indicates
//    the outside of the object.
//  --The  input polyhedra should not touch at only one point locally. If so, the current
//    implementation just ignore it (TAG SL001)
//  --Polyhedron type should be list-based or should guarantee no reallocation. We use maps
//    on pointer of halfedges,facets and vertices
//  --Polyhedral_mesh_domain requires the domain tp be closed: do not provided as input
//    an open polyhedral surface and a polyhedron with a connected component free from intersection
//OPTIMIZATIONS
//  --cdt: try using intervals? in that case, only points inside the face should be considered
//         and points on edge should be handled by hand (simply start using the point opposite to the edge)
//  --sorted_around_edge_filtered: can be done using the original supporting planes
//  --in intersection_of_Polyhedra_3: upon call to Triangle_segment_intersection_point::add_new_node, interval and exact nodes are
//    inserted into a vector. Since we do not know the final size of vector this lead to reallocation of data.
//  --in Triangle_segment_intersection_point, try using EPEC instead of Interval_nt+SC<Gmpq>
//  --use a sorted pair of indices in edge_to_hedge+simplify the code TAG_SLXX1
//  --in sew_2_marked_darts arrange how darts are passed to avoid comparing to a Point_3
//TODO:
//  --validity of the embedding: points inserted in polyhedron are approximation of the real
//    intersection points. It may happen that because of the approximation, the embedding gets
//    wrong. To avoid this, for each new triangle created, we should make an orientation test
//    with the approximated point to check if this is correct. If not, points must be moved
//    within their double interval so that all triangles incident to each of these points are correctly
//    oriented. This is probably an expensive test that can be activated only with a template parameter
//    of something similar.
//  -- We should have an option to report intersection at a point and refinement of polyhedra
//     so that the point is included into the other polyhedra. Please be careful if the point
//     already exists bool report_isolated_point or a template parameter (polyline with a unique point)
namespace CGAL
{

  namespace internal_IOP
  {

  template <class HDS, class NestedFacetConstruct, class NewNodeVertexVisitor, class PolyhedronPointPMap>
  class Triangulate_a_face : public CGAL::Modifier_base<HDS> {
    typedef typename HDS::Halfedge_handle Halfedge_handle;
    typedef typename HDS::Vertex_handle   Vertex_handle;
    typedef typename HDS::Face_handle     Face_handle;
    typedef typename HDS::Vertex          Vertex;
    typedef typename HDS::Halfedge        Halfedge;
    typedef typename HDS::Face            Face;
    typedef typename boost::property_traits<PolyhedronPointPMap>::value_type Point;

    //data members
    Face_handle current_face;
    std::map<int, Point>                                   nodes_;
    std::map<int,Vertex_handle>&                           node_to_polyhedron_vertex_;
    std::map<std::pair<int,int>,Halfedge_handle>&          edge_to_hedge_;
    std::vector<std::pair<int,int> >                       edges_to_create_;
    std::vector<CGAL::cpp11::tuple<int,int,int> >          faces_to_create_;
    NestedFacetConstruct facet_construct;
    NewNodeVertexVisitor& node_vertex_visitor;
    PolyhedronPointPMap ppmap;

    typename HDS::Halfedge::Base*
    unlock_halfedge(Halfedge_handle h){
      return static_cast<typename HDS::Halfedge::Base*>(&(*h));
    }

    typename HDS::Face::Base*
    unlock_face(Face_handle f){
      return static_cast<typename HDS::Face::Base*>(&(*f));
    }

  public:

    template <class Nodes_vector,class Triangulation>
    Triangulate_a_face( Face_handle face,
                        const Nodes_vector& nodes,
                        const std::vector<int>& node_ids,
                        std::map<int,Vertex_handle>& node_to_polyhedron_vertex,
                        std::map<std::pair<int,int>,Halfedge_handle>& edge_to_hedge,
                        const Triangulation& triangulation,
                        const NestedFacetConstruct& fc,
                        NewNodeVertexVisitor& nv,
                        PolyhedronPointPMap ppmap)
    :current_face(face), node_to_polyhedron_vertex_(node_to_polyhedron_vertex), edge_to_hedge_(edge_to_hedge), facet_construct(fc), node_vertex_visitor(nv), ppmap(ppmap)
    {
      //grab vertices to be inserted to copy them from the vector
      for (std::vector<int>::const_iterator it=node_ids.begin();it!=node_ids.end();++it)
      {
        nodes_.insert(std::make_pair(*it,nodes[*it]));
      }
      //grab edges that are not on the convex hull (these have already been created)
      for (typename Triangulation::Finite_edges_iterator
        it=triangulation.finite_edges_begin();
        it!=triangulation.finite_edges_end();
        ++it)
      {
        typename Triangulation::Vertex_handle v0=it->first->vertex((it->second+1)%3);
        typename Triangulation::Vertex_handle v1=it->first->vertex((it->second+2)%3);
        //warning in degenerate cases you can insert outsite expected convex hull edges: need exact here.
        //an alternative is to test if one the incident faces are infinite (cf assertion below)
        if ( edge_to_hedge_.find(std::make_pair(v0->info(),v1->info()))==edge_to_hedge_.end() &&
             edge_to_hedge_.find(std::make_pair(v1->info(),v0->info()))==edge_to_hedge_.end()    )
        {
          edges_to_create_.push_back( std::make_pair(v0->info(),v1->info()) );
        }
        else
            CGAL_assertion( triangulation.is_infinite(it->first->vertex(it->second)) || triangulation.is_infinite( triangulation.mirror_vertex(it->first,it->second)) );
      }
      //grab triangles.
      for (typename Triangulation::Finite_faces_iterator
        it=triangulation.finite_faces_begin();
        it!=triangulation.finite_faces_end();
        ++it)
      {
        typename Triangulation::Vertex_handle v0=it->vertex(0);
        typename Triangulation::Vertex_handle v1=it->vertex(1);
        typename Triangulation::Vertex_handle v2=it->vertex(2);
        //warning in degenerate case we can have non wanted triangles: need exact here
        faces_to_create_.push_back( CGAL::cpp11::make_tuple( v0->info(),v1->info(),v2->info() ) );
      }
    }



    void operator()( HDS& hds) {
//      std::cerr << "node_to_polyhedron_vertex_"<< std::endl;
//      for (typename std::map<int,Vertex_handle>::iterator it=node_to_polyhedron_vertex_.begin();it!=node_to_polyhedron_vertex_.end();++it)
//        std::cerr << it->first << " " << &(*(it->second)) << std::endl;

      //insert the intersection point interior to the face inside the polyhedron and
      //save their Polyhedron::vertex_handle
      for (typename std::map<int,Point>::iterator it=nodes_.begin();it!=nodes_.end();++it)
      {
        Vertex_handle v=hds.vertices_push_back(Vertex());
        node_vertex_visitor.new_vertex_added(it->first, v);
        put(ppmap, v, it->second);
        CGAL_assertion( node_to_polyhedron_vertex_.find( it->first ) == node_to_polyhedron_vertex_.end());
        node_to_polyhedron_vertex_.insert( std::make_pair(it->first,v) );
//        std::cerr << "vertices " << it->first  << " " << &(*v) << std::endl;
      }

      //insert the new halfedge and set their incident vertex
      for (typename std::vector<std::pair<int,int> >::iterator
        it=edges_to_create_.begin();it!=edges_to_create_.end();++it)
      {
        Halfedge_handle he=hds.edges_push_back(Halfedge(),Halfedge());

        //associate edge <i,j> to halfedge going from i to j with j as incident vertex
        CGAL_assertion(node_to_polyhedron_vertex_.find(it->second)!= node_to_polyhedron_vertex_.end());
        Vertex_handle v=node_to_polyhedron_vertex_.find(it->second)->second;
        unlock_halfedge(he)->set_vertex( v );
        v->set_halfedge(he);
//        std::cerr << "  --in edge " << &(*v) << std::endl;
        edge_to_hedge_.insert( std::make_pair(*it,he) );
        v=node_to_polyhedron_vertex_.find(it->first)->second;
//        std::cerr << "  --in edge " << &(*v) << std::endl;
        unlock_halfedge( he->opposite() )->set_vertex( v );
        v->set_halfedge(he->opposite());
        edge_to_hedge_.insert( std::make_pair(std::make_pair(it->second,it->first),he->opposite()) );
//        std::cerr << "edges " << it->first <<  " " << it->second << std::endl;
      }

      std::vector<CGAL::cpp11::tuple<int,int,int> >::iterator it=faces_to_create_.begin();
      Face_handle face_triangulated = current_face;
      //create the new faces and update adjacencies
      while (true)
      {
        int i=cpp11::get<0>(*it),j=cpp11::get<1>(*it),k=cpp11::get<2>(*it);
//        std::cerr << "faces " << i <<  " " << j  << " " << k<< std::endl;
        Halfedge_handle current  = edge_to_hedge_.find(std::make_pair(i,j))->second;
        Halfedge_handle next     = edge_to_hedge_.find(std::make_pair(j,k))->second;
        Halfedge_handle previous = edge_to_hedge_.find(std::make_pair(k,i))->second;


        CGAL_assertion (edge_to_hedge_.find(std::make_pair(i,j))!=edge_to_hedge_.end());
        CGAL_assertion (edge_to_hedge_.find(std::make_pair(j,k))!=edge_to_hedge_.end());
        CGAL_assertion (edge_to_hedge_.find(std::make_pair(k,i))!=edge_to_hedge_.end());

        CGAL_assertion(current->vertex()==node_to_polyhedron_vertex_.find(j)->second);
        CGAL_assertion(next->vertex()==node_to_polyhedron_vertex_.find(k)->second);
        CGAL_assertion(previous->vertex()==node_to_polyhedron_vertex_.find(i)->second);

        unlock_halfedge(current)->set_next(next);
        unlock_halfedge(next)->set_next(previous);
        unlock_halfedge(previous)->set_next(current);

        unlock_halfedge(current)->set_prev(previous);
        unlock_halfedge(next)->set_prev(current);
        unlock_halfedge(previous)->set_prev(next);

        //update face halfedge
        unlock_face(current_face)->set_halfedge(current);

        //update face of halfedges
        unlock_halfedge(current)  ->set_face(current_face);
        unlock_halfedge(next)     ->set_face(current_face);
        unlock_halfedge(previous) ->set_face(current_face);

        if ( ++it!=faces_to_create_.end() )
          current_face=hds.faces_push_back( facet_construct(*face_triangulated) );
        else
          break;
      }
    }
  };

  template <class Halfedge_handle,class Marked_set>
  Halfedge_handle
  next_marked_halfedge_around_target_vertex(Halfedge_handle h, const Marked_set& is_marked)
  {
    CGAL_assertion( is_marked.find(h)!=is_marked.end() );
    Halfedge_handle next=h->next();
    while( is_marked.find(next)==is_marked.end() )
    {
      next=next->opposite()->next();
    }
    CGAL_assertion(next!=h);
    return next;
  }

  template <class Halfedge_handle,class Marked_set>
  Halfedge_handle
  next_marked_halfedge_around_source_vertex(Halfedge_handle h, const Marked_set& is_marked)
  {
    CGAL_assertion( is_marked.find(h)!=is_marked.end() );
    Halfedge_handle prev=h->prev();
    while(is_marked.find(prev)==is_marked.end())
    {
      prev=prev->opposite()->prev();
    }
    CGAL_assertion(prev!=h);
    return prev;
  }

  } //namespace internal_IOP

template <class Polyhedron>
struct Default_facet_construct{
  typename Polyhedron::Facet operator()( const typename Polyhedron::Facet& f)
  { return f; }
};

template <class Polyhedron>
struct Default_node_vertex_visitor{
  void new_node_added(  int /* node_id */,
                        internal_IOP::Intersection_type /* type */,
                        typename Polyhedron::Halfedge_handle /* principal_edge */,
                        typename Polyhedron::Halfedge_handle /* additional_edge */,
                        bool /* is_vertex_coplanar */,
                        bool /* is_vertex_opposite_coplanar */ )
  {}

  void new_vertex_added(int /* node_id */, typename Polyhedron::Vertex_handle /* vh */){}
};

struct Default_output_builder{
  template <class T1, class T2, class T3, class T4, class T5>
  void operator()(const T1&, const T2&, const T3&, const T4&, const T5&){}
  void input_have_coplanar_facets(){}
};

template< class Polyhedron,
          class OutputBuilder_=Default,
          class Kernel_=Default,
          class EdgeMarkPropertyMap_=Default,
          class NestedFacetConstruct_=Default,
          class NewNodeVertexVisitor_=Default,
          class PolyhedronPointPMap_=Default
        >
class Node_visitor_refine_polyhedra{
//Default typedefs
  typedef typename Default::Get<OutputBuilder_, Default_output_builder >::type OutputBuilder;
  typedef typename Default::Get<EdgeMarkPropertyMap_, Corefinement::Dummy_edge_mark_property_map<Polyhedron> >::type EdgeMarkPropertyMap;
  typedef typename Default::Get<NestedFacetConstruct_, Default_facet_construct<Polyhedron > >::type NestedFacetConstruct;
  typedef typename Default::Get<NewNodeVertexVisitor_, Default_node_vertex_visitor<Polyhedron> >::type NewNodeVertexVisitor;
  typedef typename Default::Get<PolyhedronPointPMap_, Default_polyhedron_ppmap<Polyhedron> >::type PolyhedronPointPMap;
  typedef typename Default::Get<Kernel_, typename Kernel_traits< typename boost::property_traits<PolyhedronPointPMap>::value_type >::Kernel >::type Kernel;
//typedefs
  typedef typename Polyhedron::Halfedge_handle                         Halfedge_handle;
  typedef typename Polyhedron::Halfedge_const_handle                   Halfedge_const_handle;
  typedef typename Polyhedron::Face_handle                             Face_handle;
  typedef typename Polyhedron::Halfedge                                Halfedge;
  typedef typename Polyhedron::Vertex_handle                           Vertex_handle;
  typedef internal_IOP::Compare_handles<Polyhedron,CGAL::Tag_false>    Cmp_handle; //This ensures uniqueness of edges when comparing halfedges
  typedef internal_IOP::Compare_unik_address<Polyhedron>               Cmp_unik_ad; //This ensures uniqueness of edges when comparing halfedges

  //constrained triangulation used for triangulation interior of faces
  #ifdef DO_NO_USE_EXACT_CDT
  typedef CGAL::Triangulation_vertex_base_with_info_2<int,Kernel>       Vbi;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi,Fb>                  TDS_2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,CGAL::No_intersection_tag> CDT;  //DO WE NEED DELAUNAY????
  #else
  /// \todo change this, use it only if not already exact
  typedef CGAL::Exact_predicates_exact_constructions_kernel             Exact_kernel;
  typedef CGAL::Triangulation_vertex_base_with_info_2<int,Exact_kernel> Vbi;
  typedef CGAL::Constrained_triangulation_face_base_2<Exact_kernel>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi,Fb>                  TDS_2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Exact_kernel,TDS_2,CGAL::No_intersection_tag> CDT;  //DO WE NEED DELAUNAY????
  #endif

  typedef std::map<Halfedge_handle,Polyhedron*,Cmp_unik_ad>            Hedge_to_polyhedron_map;

  typedef std::vector<int>                                             Node_ids;
  typedef std::map< Face_handle,Node_ids,Cmp_handle >                  In_face_map;
  typedef std::map< Halfedge_handle,Node_ids,Cmp_unik_ad >             In_halfedge_map;
  //to keep the correspondance between node_id and vertex_handle in each polyhedron
  typedef std::map<int,Vertex_handle> Node_to_polyhedron_vertex_map;
  typedef std::map<Polyhedron*, Node_to_polyhedron_vertex_map > Poly_to_map_node;
  //to maintain an polyhedron halfedge on each polyline + pair<bool,int>
  //with first = "is the key (pair<int,int>) was reversed?" and second is the number of edges +1 in the polyline
  typedef std::map< std::pair<int,int>, std::pair< std::map<Polyhedron*,Halfedge_handle>,std::pair<bool,int> > > An_edge_per_polyline_map;
  //to handle coplanar halfedge of polyhedra that are full in the intersection
  typedef std::map< int,Halfedge_handle >                              Node_to_target_of_hedge_map;
  typedef std::map< Polyhedron*,Node_to_target_of_hedge_map>           Poly_to_vertices_on_intersection_map;

//data members
  Hedge_to_polyhedron_map               hedge_to_polyhedron;
  In_face_map                           in_face;
  In_halfedge_map                       in_hedge;
  std::map< int,std::set<int> >         graph_of_constraints;
  boost::dynamic_bitset<>               is_node_of_degree_one;
  std::map< int,std::set<int> >         coplanar_constraints;
  An_edge_per_polyline_map              an_edge_per_polyline;
  typename An_edge_per_polyline_map::iterator last_polyline;
  Poly_to_vertices_on_intersection_map  poly_to_vertices_on_inter;
  Poly_to_map_node                      polyhedron_to_map_node_to_polyhedron_vertex;
  std::set<int>                         non_manifold_nodes; //contain nodes that are original vertices of input polyhedron and that neighborhood is not a topological disk
  std::map<Vertex_handle,int>           nodes_that_are_original_vertices;//to keep the correspondance between original polyhedron vertices that are also nodes

  OutputBuilder                         output_builder;

  //   new_hedge    hedge
  //  ----------->   ----------->
  //               v
  //  <-----------   <-----------
  //   new_opposite     opposite
  //
  template <class Nodes_vector>
  Vertex_handle split_edge( Halfedge_handle hedge,
                   int node_id,
                   const Nodes_vector& nodes,
                   Polyhedron& P)
  {
    internal_IOP::Split_halfedge<typename Polyhedron::HalfedgeDS> delegated(hedge);
    P.delegate( delegated );
    CGAL_assertion(P.is_valid());

    Vertex_handle vh=boost::prior(P.vertices_end());
    node_vertex_visitor.new_vertex_added(node_id, vh);
    put(ppmap, vh, nodes[node_id]);
    CGAL_assertion(get(ppmap,vh)==nodes[node_id]);

    //update marker tags. If the edge was marked, then the resulting edges in the split must be marked
    if ( get(m_edge_mark_pmap,std::make_pair(hedge,&P)) )
    {
      CGAL_assertion( get(m_edge_mark_pmap,std::make_pair(hedge->opposite(),&P)) );
      put(m_edge_mark_pmap,std::make_pair(hedge->prev(),&P),true);
      put(m_edge_mark_pmap,std::make_pair(hedge->prev()->opposite(),&P),true);
    }

    return vh;
  }

  //sort node ids so that we can split the hedge
  //consecutively
  template <class Nodes_vector>
  void sort_vertices_along_hedge(std::vector<int>& node_ids,Halfedge_handle hedge,const Nodes_vector& nodes)
  {
    std::sort(node_ids.begin(),
              node_ids.end(),
              internal_IOP::Order_along_a_halfedge<Polyhedron,PolyhedronPointPMap, Nodes_vector,Is_polyhedron_const>(hedge,nodes, ppmap)
    );
  }

  //insert intersection as constrained edges  in a CDT triangulation
  template <class CDT>
  void insert_constrained_edges_coplanar_case(int node_id,
                                              CDT& triangulation,
                                              std::map<int,typename CDT::Vertex_handle>& id_to_CDT_vh)
  {
    if (node_id < number_coplanar_vertices){
      //XSL_TAG_CPL_VERT
      //Insert constrained edges from coplanar facets that have been retriangulated. This ensure that triangulations are compatible
      std::map< int,std::set<int> >::iterator it_neighbors=coplanar_constraints.find(node_id);
      if (it_neighbors!=coplanar_constraints.end())
      {
        typename CDT::Vertex_handle vh=id_to_CDT_vh.find(node_id)->second;
        for (std::set<int>::iterator it_n=it_neighbors->second.begin();it_n!=it_neighbors->second.end();++it_n){
          typename std::map<int,typename CDT::Vertex_handle>::iterator it_vh=id_to_CDT_vh.find(*it_n);
          // this condition ensures to consider only graph edges that are in the same triangle (not in a neighbor one when involving node on a triangle edge)
          // here we can't make the difference between a point on the interior or the boundary, so points_on_triangle is not used.
          if ( it_vh!=id_to_CDT_vh.end() ){
            triangulation.insert_constraint(vh,id_to_CDT_vh.find(*it_n)->second);
          }
        }
      }
    }
  }
  //insert intersection as constrained edges  in a CDT triangulation
  template <class CDT,class Constrained_edges_map>
  void insert_constrained_edges(Node_ids& node_ids, //index of vertices we are interested in
                                CDT& triangulation,
                                std::map<int,typename CDT::Vertex_handle>& id_to_CDT_vh,
                                Constrained_edges_map& constrained_edges, //list of pair of int to indicate edges that are constrained
                                bool points_on_triangle=false)
  {
    for (Node_ids::iterator it_node_id=node_ids.begin();it_node_id!=node_ids.end();++it_node_id){
      std::map< int,std::set<int> >::iterator it_neighbors=graph_of_constraints.find(*it_node_id);
      if (it_neighbors!=graph_of_constraints.end())
      {
        typename CDT::Vertex_handle vh=id_to_CDT_vh.find(*it_node_id)->second;
        for (std::set<int>::iterator it_n=it_neighbors->second.begin();it_n!=it_neighbors->second.end();++it_n){
          typename std::map<int,typename CDT::Vertex_handle>::iterator it_vh=id_to_CDT_vh.find(*it_n);
          // this condition ensures to consider only graph edges that are in the same triangle (not in a neighbor one when involving node on a triangle edge)
          if ( !points_on_triangle || it_vh!=id_to_CDT_vh.end() ){
            CGAL_assertion(it_vh!=id_to_CDT_vh.end());
            triangulation.insert_constraint(vh,id_to_CDT_vh.find(*it_n)->second);
            constrained_edges.push_back(std::make_pair(*it_node_id,*it_n));
          }
        }
      }
      #ifdef CGAL_COREFINEMENT_DEBUG
      else
      {
        std::cout << "X0: Found an isolated point" << std::endl;
      }
      #endif

      insert_constrained_edges_coplanar_case(*it_node_id,triangulation,id_to_CDT_vh);
    }
  }

  std::pair<int,int> make_sorted_pair(int i,int j) const {return i<j ? std::make_pair(i,j):std::make_pair(j,i);}

  void update_edge_per_polyline(Polyhedron* P,std::pair<int,int> indices,typename Polyhedron::Halfedge_handle hedge)
  {
    std::pair<int,int> sorted_pair=make_sorted_pair(indices.first,indices.second);
    typename An_edge_per_polyline_map::iterator it=an_edge_per_polyline.find(sorted_pair);
    if (it!=an_edge_per_polyline.end()){
      it->second.first.insert(std::make_pair( P,sorted_pair.first==indices.first?hedge:hedge->opposite() ));
    }
  }

  //keep track of the fact that a polyhedron original vertex is a node
  void all_incident_faces_got_a_node_as_vertex(Halfedge_handle incident_to_vertex_edge,int node_id)
  {
    nodes_that_are_original_vertices.insert(std::make_pair(incident_to_vertex_edge->vertex(),node_id));
  }

  //if an original polyhedron vertex is also a node, do no use a fake id
  void set_triangle_boundary_indices(
    Vertex_handle* triangle_boundary,
    int* triangle_boundary_indices)
  {
    triangle_boundary_indices[0]=-1;
    triangle_boundary_indices[1]=-2;
    triangle_boundary_indices[2]=-3;

    for (int k=0;k<3;++k){
      typename std::map<Vertex_handle,int>::iterator it=nodes_that_are_original_vertices.find(triangle_boundary[k]);
      if (it!=nodes_that_are_original_vertices.end())
        triangle_boundary_indices[k]=it->second;
    }
  }

  int number_coplanar_vertices; //number of intersection points between coplanar facets, see fixes XSL_TAG_CPL_VERT
  EdgeMarkPropertyMap m_edge_mark_pmap;     //property map to mark halfedge of the original polyhedra that are on the intersection
  NestedFacetConstruct facet_construct;  // functor called to create new triangular faces inside a given face
  NewNodeVertexVisitor node_vertex_visitor; // functor called when a new node is created and when a new vertex is added
  PolyhedronPointPMap ppmap;
public:
  Node_visitor_refine_polyhedra (
    OutputBuilder output_builder_=OutputBuilder(),
    PolyhedronPointPMap ppmap = PolyhedronPointPMap(),
    EdgeMarkPropertyMap pmap=EdgeMarkPropertyMap(),
    const NestedFacetConstruct& fc = NestedFacetConstruct(),
    const NewNodeVertexVisitor& nv = NewNodeVertexVisitor() )
  : output_builder(output_builder_)
  , m_edge_mark_pmap(pmap)
  , facet_construct(fc)
  , node_vertex_visitor(nv)
  , ppmap(ppmap)
  {}

  typedef internal_IOP::Predicates_on_constructions Node_storage_type;
  typedef Tag_false Is_polyhedron_const;
  static const bool do_need_vertex_graph = true;  //because we need to know which edges are constrained

  void set_number_of_intersection_points_from_coplanar_facets(int n){
    number_coplanar_vertices=n;
  }

  void input_have_coplanar_facets()
  {
    output_builder.input_have_coplanar_facets();
  }

  void check_node_on_non_manifold_vertex(int node_id,Halfedge_handle hedge){
    //we turn around the hedge and check no halfedge is a border halfedge
    Halfedge_handle curr=hedge;
    do{
      if ( curr->is_border_edge() ){
        non_manifold_nodes.insert(node_id);
        return;
      }
      curr=curr->next()->opposite();
    }
    while(curr!=hedge);
  }

  void check_node_on_non_manifold_edge(int node_id,Halfedge_handle hedge){
    if ( hedge->is_border_edge() ) non_manifold_nodes.insert(node_id);
  }

  void new_node_added(int node_id,
                      internal_IOP::Intersection_type type,
                      Halfedge_handle principal_edge,
                      Halfedge_handle additional_edge,
                      bool is_vertex_coplanar,
                      bool is_vertex_opposite_coplanar)
  {
    //forward to the visitor
    node_vertex_visitor.new_node_added(node_id, type, principal_edge, additional_edge, is_vertex_coplanar, is_vertex_opposite_coplanar);
    switch(type)
    {
      case internal_IOP::FACET: //Facet intersected by an edge
      {
        typename In_face_map::iterator it_fmap=in_face.insert(std::make_pair(additional_edge->face(), Node_ids())).first;
        it_fmap->second.push_back(node_id);
      }
      break;
      case internal_IOP::EDGE: //Edge intersected by an edge
      {
        typename In_halfedge_map::iterator it_hedge_map=in_hedge.insert(std::make_pair(additional_edge,Node_ids())).first;
        it_hedge_map->second.push_back(node_id);
        check_node_on_non_manifold_edge(node_id,additional_edge);
      }
      break;
      case internal_IOP::VERTEX:
      {
        //grab original vertex that is on commom intersection
        typename Hedge_to_polyhedron_map::iterator it=hedge_to_polyhedron.find(additional_edge->facet()->halfedge());
        CGAL_assertion(it!=hedge_to_polyhedron.end());
        poly_to_vertices_on_inter[it->second].insert(std::make_pair(node_id,additional_edge));
        polyhedron_to_map_node_to_polyhedron_vertex[it->second].insert(std::make_pair(node_id,additional_edge->vertex()));
        all_incident_faces_got_a_node_as_vertex(additional_edge,node_id);
        check_node_on_non_manifold_vertex(node_id,additional_edge);
      }
      break;
      default:
      return;
    }

    CGAL_assertion(!is_vertex_coplanar || !is_vertex_opposite_coplanar); //coplanar edge are not forwarded


    if ( is_vertex_coplanar )
    {
      //grab original vertex that is on commom intersection
      typename Hedge_to_polyhedron_map::iterator it=hedge_to_polyhedron.find(principal_edge->facet()->halfedge());
      poly_to_vertices_on_inter[it->second].insert(std::make_pair(node_id,principal_edge));
      polyhedron_to_map_node_to_polyhedron_vertex[it->second].insert(std::make_pair(node_id,principal_edge->vertex()));
      all_incident_faces_got_a_node_as_vertex(principal_edge,node_id);
      check_node_on_non_manifold_vertex(node_id,principal_edge);
    }
    else{
      if ( is_vertex_opposite_coplanar ){
        //grab original vertex that is on commom intersection
        typename Hedge_to_polyhedron_map::iterator it=hedge_to_polyhedron.find(principal_edge->facet()->halfedge());
        poly_to_vertices_on_inter[it->second].insert(std::make_pair(node_id,principal_edge->opposite()));
        polyhedron_to_map_node_to_polyhedron_vertex[it->second].insert(std::make_pair(node_id,principal_edge->opposite()->vertex()));
        all_incident_faces_got_a_node_as_vertex(principal_edge->opposite(),node_id);
        check_node_on_non_manifold_vertex(node_id,principal_edge->opposite());
      }
      else{
        //handle intersection on principal edge
        typename In_halfedge_map::iterator it_hedge_map=in_hedge.insert(std::make_pair(principal_edge,Node_ids())).first;
        it_hedge_map->second.push_back(node_id);
        check_node_on_non_manifold_edge(node_id,principal_edge);
      }
    }
  }

  template<class Iterator>
  void annotate_graph(Iterator begin,Iterator end)
  {
//    std::cout << "Annotation graph..." << std::endl;
    int node_id = 0;
    is_node_of_degree_one.resize(std::distance(begin, end));
    for (Iterator it=begin;it!=end;++it, ++node_id)
    {
      if (non_manifold_nodes.count(node_id)) it->make_terminal();
      const std::set<int>& neighbors = it->neighbors;
      graph_of_constraints.insert(std::make_pair(node_id,neighbors));
      if (neighbors.size()==1)
        is_node_of_degree_one.set(node_id);
    }
  }

  void update_terminal_nodes(std::vector<bool>&)
  {
    CGAL_assertion(!"Must not call this function");
  }

  void add_filtered_intersection(Halfedge_handle eh,Halfedge_handle fh,Polyhedron& Pe,Polyhedron& Pf){
    //use the representant halfedge of the facet as key
    //--set polyhedron for the two facets incident to the edge
    CGAL_assertion(!eh->is_border());
    hedge_to_polyhedron.insert(std::make_pair(eh->facet()->halfedge(),&Pe));
    if ( !eh->opposite()->is_border() )
      hedge_to_polyhedron.insert(std::make_pair(eh->opposite()->facet()->halfedge(),&Pe));
    //--set polyhedron for the facet intersected by the edge
    hedge_to_polyhedron.insert(std::make_pair(fh->facet()->halfedge(),&Pf));
  }


  struct Polyhedron_face_boundary{
    std::vector<int> node_ids_array[3]; // the node_ids on each halfedges
    std::map<Halfedge_handle,int,Cmp_unik_ad> hedges_ids;
    Halfedge_handle halfedges[3]; //the three halfedges of the original face
    Vertex_handle   vertices[3];  //the three vertices  of the original face
    //node_ids_array[0] corresponds to the original edge vertices[0],vertices[1] = halfedges[0]
    //node_ids_array[1] corresponds to the original edge vertices[1],vertices[2] = halfedges[1]
    //node_ids_array[2] corresponds to the original edge vertices[2],vertices[0] = halfedges[2]
    Polyhedron_face_boundary(Halfedge_handle first)
    {
      CGAL_assertion(first->next()->next()->next()==first); //the face is a triangle
      hedges_ids.insert(std::make_pair(first,0));
      hedges_ids.insert(std::make_pair(first->next(),1));
      hedges_ids.insert(std::make_pair(first->next()->next(),2));
      halfedges[0]=first;
      halfedges[1]=first->next();
      halfedges[2]=first->next()->next();

      vertices[0]=halfedges[0]->opposite()->vertex();
      vertices[1]=halfedges[1]->opposite()->vertex();
      vertices[2]=halfedges[2]->opposite()->vertex();
    }

    //used when object was created with hedge but opposite was used to split the original face
    void update_original_halfedge(Halfedge_handle original,Halfedge_handle new_hedge)
    {
      typename std::map<Halfedge_handle,int,Cmp_unik_ad>::iterator it_id=hedges_ids.find(original);
      CGAL_assertion(it_id!=hedges_ids.end());
      int index=it_id->second;
      CGAL_assertion(halfedges[index]==original);
      hedges_ids.erase(it_id);
      hedges_ids.insert(std::make_pair(new_hedge,index));
      halfedges[index]=new_hedge;
    }

    template <class Iterator>
    void copy_node_ids(Halfedge_handle hedge,Iterator begin,Iterator end)
    {
      typename std::map<Halfedge_handle,int,Cmp_unik_ad>::iterator it_id=hedges_ids.find(hedge);
      CGAL_assertion(it_id!=hedges_ids.end());
      std::copy(begin,end,std::back_inserter(node_ids_array[it_id->second]));
    }
  };


  void start_new_polyline(int i, int j)
  {
    if ( i==j ) //case of a single point
    {
      //TAG SL001
      //nothing is done
      return;
    }
    std::pair<typename An_edge_per_polyline_map::iterator,bool> res=
      an_edge_per_polyline.insert(
        std::make_pair( make_sorted_pair(i,j),
          std::make_pair( std::map<Polyhedron*,Halfedge_handle>(),std::make_pair(false,0))  )
      );
    CGAL_assertion(res.second);
    last_polyline=res.first;
    if ( i !=last_polyline->first.first )
      last_polyline->second.second.first=true;
  }

  void add_node_to_polyline(int){
    ++(last_polyline->second.second.second);
  }

  void new_input_polyhedron(Polyhedron& P)
  {
    typedef std::pair<typename Poly_to_map_node::iterator,bool> Res;
    CGAL_USE_TYPE(Res);
    CGAL_assertion_code(Res res = )
      polyhedron_to_map_node_to_polyhedron_vertex.insert(std::make_pair( &P,Node_to_polyhedron_vertex_map() ));
    CGAL_assertion(res.second == true);
  }

  //1) split_halfedges and retriangulate faces with no intersection point interior to the facet
  //2) retriangulate using a constrained Delaunay triangulation each triangle in each Polyhedron that contains at least
  //   one intersection point inside the facet
  //3) mark polyhedron edges that are on the intersection
  //4) create one output polyhedron per connected component of polyhedron, connected by an edge which is not an intersection edge
  //5) import each piece into a common combinatorial map
  //6) glue all the pieces together
  template <class Nodes_vector>
  void finalize(const Nodes_vector& nodes){
    //mark halfedge that are on the intersection
    //SL: I needed to use a map because to get the orientation around the edge,
    //    I need to know in the case the third vertex is a node its index (for exact construction)
    typedef std::map<Halfedge_const_handle,std::pair<int,int>,Cmp_unik_ad > Border_halfedges_map;
    Border_halfedges_map border_halfedges;

    //store for each triangle facet which boundary is intersected by the other surface,
    //original vertices (and halfedges in the refined mesh pointing on these vertices)
    typedef std::map<Face_handle,Polyhedron_face_boundary,Cmp_handle> Faces_boundary;
    Faces_boundary faces_boundary;

    //0) For each polyhedron, collect original vertices that belongs to the intersection.
    //   From the graph of constraints, extract intersection edges that are incident to such vertices. In case
    //   there exists another original vertex adjacent to the first one found, this halfedge must be
    //   marked on the boundary (and possibly update an_edge_per_polyline).
    //   This is done first to avoid halfedges stored to be modified in the steps following.
    for (typename Poly_to_vertices_on_intersection_map::iterator
      it=poly_to_vertices_on_inter.begin();
      it!=poly_to_vertices_on_inter.end();
      ++it)
    {
      Polyhedron* poly=it->first;
      std::set<std::pair<int,int> > already_done;
      Node_to_target_of_hedge_map& nodes_to_hedge=it->second;
      for(typename Node_to_target_of_hedge_map::iterator
        it_node_2_hedge=nodes_to_hedge.begin();
        it_node_2_hedge!=nodes_to_hedge.end();
        ++it_node_2_hedge)
      {
        int node_id_of_first=it_node_2_hedge->first;
        std::map< int,std::set<int> >::iterator it_neighbors=graph_of_constraints.find(node_id_of_first);
        if ( it_neighbors!=graph_of_constraints.end() )
        {
          std::set<int>& neighbors=it_neighbors->second;
          for (std::set<int>::iterator it_id=neighbors.begin();it_id!=neighbors.end();++it_id){
            if ( already_done.find(std::make_pair(*it_id,node_id_of_first))!=already_done.end() ) continue;//already done for the opposite
            typename Node_to_target_of_hedge_map::iterator it_node_2_hedge_two=nodes_to_hedge.find(*it_id);
            if ( it_node_2_hedge_two!=nodes_to_hedge.end() ) //a full edge is on intersection
            {
              //get the corresponding halfedge with vertex corresponding to node_id_of_first
              Halfedge_handle hedge=it_node_2_hedge->second;
              CGAL_assertion_code(Halfedge_handle start=hedge;)
              while ( hedge->opposite()->vertex()!=it_node_2_hedge_two->second->vertex() ){
                hedge=hedge->next()->opposite();
                CGAL_assertion(hedge!=start);
              }
              std::pair<int,int> edge_pair(*it_id,node_id_of_first);
              if ( border_halfedges.insert( std::make_pair(hedge,edge_pair) ).second)
              {
                put(m_edge_mark_pmap,std::make_pair(hedge,poly),true);
                put(m_edge_mark_pmap,std::make_pair(hedge->opposite(),poly),true);
              }
              update_edge_per_polyline(poly,edge_pair,hedge);
              //save the fact that we already handle this edge
              already_done.insert(std::make_pair(node_id_of_first,*it_id));
            }
          }
        }
        #ifdef CGAL_COREFINEMENT_DEBUG
        else
        {
          std::cout << "X1: Found an isolated point" << std::endl;
        }
        #endif
      }
    }

    //1) First split halfedges cut by the intersection polyline(s)
    for (typename In_halfedge_map::iterator it=in_hedge.begin();it!=in_hedge.end();++it)
    {
      Halfedge_handle hedge=it->first; //the halfedge to be split (and its opposite too)
      Node_ids& node_ids = it->second;   //indices of the intersection points to be inserted
      CGAL_assertion( std::set<int>(node_ids.begin(), node_ids.end()).size()==node_ids.size() );
      typename Hedge_to_polyhedron_map::iterator  it_poly=hedge_to_polyhedron.find( hedge->facet()->halfedge() );
      CGAL_assertion(it_poly!=hedge_to_polyhedron.end());
      Polyhedron* P=it_poly->second;  //the polyhedron in which vertices should be added

      sort_vertices_along_hedge(node_ids,hedge,nodes);

      //save original face and nodes for face of hedge (1)
      if ( !hedge->is_border() ){
        typename Faces_boundary::iterator it_face=faces_boundary.find(hedge->face());
        if (it_face==faces_boundary.end())
          it_face=faces_boundary.insert(std::make_pair(hedge->face(),Polyhedron_face_boundary(hedge))).first;
        it_face->second.copy_node_ids(hedge,node_ids.begin(),node_ids.end());
      }

      //save original face and nodes for face of hedge->opposite (2)
      typename Faces_boundary::iterator opposite_original_info=faces_boundary.end();
      if ( !hedge->opposite()->is_border() ){
        opposite_original_info=faces_boundary.find(hedge->opposite()->face());
        if (opposite_original_info==faces_boundary.end())
          opposite_original_info=faces_boundary.insert(std::make_pair(hedge->opposite()->face(),Polyhedron_face_boundary(hedge->opposite()))).first;
        opposite_original_info->second.copy_node_ids(hedge->opposite(),node_ids.rbegin(),node_ids.rend());
      }

      typename Poly_to_map_node::iterator it_map=polyhedron_to_map_node_to_polyhedron_vertex.find(P);
      CGAL_assertion(it_map!=polyhedron_to_map_node_to_polyhedron_vertex.end());
      //a map to identify the vertex in the polyhedron corresponding to an intersection point
      Node_to_polyhedron_vertex_map& node_to_polyhedron_vertex=it_map->second;

      CGAL_assertion_code(Vertex_handle original_vertex=hedge->opposite()->vertex();)

      //We need an edge incident to the source vertex of hedge. This is the first opposite edge created.
      bool first=true; Halfedge_handle hedge_incident_to_src;
      //do split the edges
      for (std::vector<int>::const_iterator it_id=node_ids.begin();it_id!=node_ids.end();++it_id){
        Vertex_handle v=split_edge(hedge, *it_id, nodes, *P);
        node_to_polyhedron_vertex.insert(std::make_pair(*it_id,v));
        if (first){
          first=false;
          hedge_incident_to_src=hedge->opposite()->next();
        }
      }

      CGAL_assertion(hedge_incident_to_src->vertex()==original_vertex);
      CGAL_assertion(hedge_incident_to_src->face()==hedge->opposite()->face());

      //save original face and nodes for face of hedge->opposite (2)
      if ( !hedge->opposite()->is_border() ){
        CGAL_assertion(opposite_original_info!=faces_boundary.end());
        opposite_original_info->second.update_original_halfedge(hedge->opposite(),hedge_incident_to_src);
      }

      //insert the two incident faces in in_face map so that they will be triangulated.
      if (!hedge->is_border()) in_face.insert(std::make_pair(hedge->face(),Node_ids()));
      if (!hedge->opposite()->is_border()) in_face.insert(std::make_pair(hedge->opposite()->face(),Node_ids()));
    }

    //2)triangulation of the triangle faces containing intersection point in their interior
    //  and also those with intersection points only on the boundary.
    for (typename In_face_map::iterator it=in_face.begin();it!=in_face.end();++it)
    {
      Face_handle f = it->first; //the face to be retriangulated
      Node_ids& node_ids  = it->second; //the index of the intersection point that are interior to the face
      CGAL_assertion(std::set<int>(node_ids.begin(), node_ids.end()).size()==node_ids.size());
      typename Faces_boundary::iterator it_fb=faces_boundary.find(f);


      typename Hedge_to_polyhedron_map::iterator it_polyhedron = hedge_to_polyhedron.find (f->halfedge()); //we can do this because the halfedge is still the same (at least its address)+no Face::set_halfedge called
      CGAL_assertion(it_polyhedron != hedge_to_polyhedron.end());
      Polyhedron* P=it_polyhedron->second;
      typename Poly_to_map_node::iterator it_map=polyhedron_to_map_node_to_polyhedron_vertex.find(P);
      CGAL_assertion(it_map!=polyhedron_to_map_node_to_polyhedron_vertex.end());
      //a map to identify the vertex in the polyhedron corresponding to an intersection point
      Node_to_polyhedron_vertex_map& node_to_polyhedron_vertex=it_map->second;

      std::map<int,typename CDT::Vertex_handle> id_to_CDT_vh;

      //associate an edge of the triangulation to a halfedge in a given polyhedron
      std::map<std::pair<int,int>,Halfedge_handle> edge_to_hedge;

      Vertex_handle triangle_boundary[3];
      int triangle_boundary_indices[3]; //the node_id of the triangle original vertex or a fake id
      if (it_fb!=faces_boundary.end()){ //the boundary of the triangle face was refined
        triangle_boundary[0]=it_fb->second.vertices[0];
        triangle_boundary[1]=it_fb->second.vertices[1];
        triangle_boundary[2]=it_fb->second.vertices[2];
        set_triangle_boundary_indices(triangle_boundary,triangle_boundary_indices);
      }
      else{
        triangle_boundary[0]=f->halfedge()->vertex(); //-1
        triangle_boundary[1]=f->halfedge()->next()->vertex(); //-2
        triangle_boundary[2]=f->halfedge()->next()->next()->vertex(); //-3
        CGAL_assertion(f->halfedge()->next()->next()->next()==f->halfedge());//check this is a triangle
        set_triangle_boundary_indices(triangle_boundary,triangle_boundary_indices);
        edge_to_hedge.insert (std::make_pair( std::make_pair( triangle_boundary_indices[2],triangle_boundary_indices[0] ) , f->halfedge() ) );
        edge_to_hedge.insert (std::make_pair( std::make_pair( triangle_boundary_indices[0],triangle_boundary_indices[1] ) , f->halfedge()->next() ) );
        edge_to_hedge.insert (std::make_pair( std::make_pair( triangle_boundary_indices[1],triangle_boundary_indices[2] ) , f->halfedge()->next()->next() ) );
      }


      #ifdef DO_NO_USE_EXACT_CDT
      typename Kernel::Plane_3 plane( get(ppmap,triangle_boundary[0]),get(ppmap,triangle_boundary[1]),get(ppmap,triangle_boundary[2]));
      #else
      CGAL::Cartesian_converter<Kernel,Exact_kernel> convert;
      typename Exact_kernel::Plane_3 plane(convert(get(ppmap,triangle_boundary[0])),convert(get(ppmap,triangle_boundary[1])),convert(get(ppmap,triangle_boundary[2])));
      #endif
      CDT triangulation;
      //insert point inside face
      for (std::vector<int>::iterator it_node_id=node_ids.begin();it_node_id!=node_ids.end();++it_node_id){
        #ifdef DO_NO_USE_EXACT_CDT
        typename CDT::Vertex_handle vh=triangulation.insert(plane.to_2d(nodes[*it_node_id]));
        #else
        typename CDT::Vertex_handle vh=triangulation.insert(plane.to_2d(nodes.exact_node(*it_node_id)));
        #endif
        vh->info()=*it_node_id;
        id_to_CDT_vh.insert(std::make_pair(*it_node_id,vh));
      }


      typename CDT::Vertex_handle triangle_vertices[3];
      #ifdef DO_NO_USE_EXACT_CDT
      triangle_vertices[0]=triangulation.insert(plane.to_2d(get(ppmap,triangle_boundary[0])));
      triangle_vertices[1]=triangulation.insert(plane.to_2d(get(ppmap,triangle_boundary[1])));
      triangle_vertices[2]=triangulation.insert(plane.to_2d(get(ppmap,triangle_boundary[2])));
      #else
      //we can do this because these are input points.
      triangle_vertices[0]=triangulation.insert(plane.to_2d(convert(get(ppmap,triangle_boundary[0]))));
      triangle_vertices[1]=triangulation.insert(plane.to_2d(convert(get(ppmap,triangle_boundary[1]))));
      triangle_vertices[2]=triangulation.insert(plane.to_2d(convert(get(ppmap,triangle_boundary[2]))));
      #endif

      triangle_vertices[0]->info()=triangle_boundary_indices[0];
      triangle_vertices[1]->info()=triangle_boundary_indices[1];
      triangle_vertices[2]->info()=triangle_boundary_indices[2];
      //insert face_extremities: we use operator[] because indice -1,-2,-3 are used in each loop and are specific to the current face
      node_to_polyhedron_vertex[-1]=triangle_boundary[0];
      node_to_polyhedron_vertex[-2]=triangle_boundary[1];
      node_to_polyhedron_vertex[-3]=triangle_boundary[2];

      //if one of the triangle original vertex is also a node
      for (int ik=0;ik<3;++ik){
        if ( triangle_boundary_indices[ik]>=0 )
          id_to_CDT_vh.insert(std::make_pair(triangle_boundary_indices[ik],triangle_vertices[ik]));
      }
      //insert points on edges
      #ifdef DO_NO_USE_EXACT_CDT
      //and constrains these edges
      #endif
      if (it_fb!=faces_boundary.end()) //is there at least one intersection point on the boundary of the face?
      {
        //in the following loop, for each original edge of the triangle, we insert the constrained edges
        // and we recover the halfedge_handle corresponding to these constrained (they are already in the polyhedron)
        for (int i=0;i<3;++i){
//          std::cerr << "Boundary edges" << std::endl;
//          std::cerr <<  "  " << -1-i <<std::endl;
          //handle case of halfedge starting at triangle_vertices[i] and ending at triangle_vertices[(i+1)%3]
          Node_ids& bounding_ids=it_fb->second.node_ids_array[i];
          typename CDT::Vertex_handle previous=triangle_vertices[i];
          int previous_index=triangle_boundary_indices[i]; //index of original Polyhedron vertex
          Halfedge_handle hedge = it_fb->second.halfedges[ (i+2) % 3]->next();
          CGAL_assertion( hedge->opposite()->vertex()==it_fb->second.vertices[i] );
          if (!bounding_ids.empty()){ //is there al least one intersection point on this edge?
            for (Node_ids::iterator it_id=bounding_ids.begin();it_id!=bounding_ids.end();++it_id){
//              std::cerr << "  "<<  *it_id << std::endl;
              #ifdef DO_NO_USE_EXACT_CDT
              typename CDT::Vertex_handle vh=triangulation.insert(plane.to_2d(nodes[*it_id]));
              #else
              typename CDT::Vertex_handle vh=triangulation.insert(plane.to_2d(nodes.exact_node(*it_id)));
              #endif
              vh->info()=*it_id;
              id_to_CDT_vh.insert(std::make_pair(*it_id,vh));
              #ifdef DO_NO_USE_EXACT_CDT
              triangulation.insert_constraint(previous,vh);
              #endif
              edge_to_hedge.insert (std::make_pair( std::make_pair(previous_index,*it_id),hedge) );
              previous=vh;
              hedge=hedge->next();
              previous_index=*it_id;
            }
          }
          else{
            CGAL_assertion( it_fb->second.halfedges[i]->vertex() == it_fb->second.vertices[ (i+1) % 3 ] );
            CGAL_assertion( it_fb->second.halfedges[i]->opposite()->vertex() == it_fb->second.vertices[ i ] );
          }
          CGAL_assertion(hedge==it_fb->second.halfedges[i]);
          edge_to_hedge.insert (std::make_pair( std::make_pair(previous_index,triangle_boundary_indices[(i+1) % 3]) , it_fb->second.halfedges[i] ) );
//          std::cerr <<  "  " << -1 - ( (i+1) % 3 ) <<std::endl;
          #ifdef DO_NO_USE_EXACT_CDT
          triangulation.insert_constraint(previous,triangle_vertices[(i+1)%3]);
          #endif
        }
      }

      std::list<std::pair<int,int> > constrained_edges;

      //insert constraints that are interior to the triangle (in the case no edges are collinear in the meshes)
      insert_constrained_edges(node_ids,triangulation,id_to_CDT_vh,constrained_edges);

      //insert constraints between points that are on the boundary (not a contrained on the triangle boundary)
      if (it_fb!=faces_boundary.end()) //is there at least one intersection point on the boundary of the face?
      {
        for (int i=0;i<3;++i){
          Node_ids& bounding_ids=it_fb->second.node_ids_array[i];
          insert_constrained_edges(bounding_ids,triangulation,id_to_CDT_vh,constrained_edges,true);
        }
      }

      //insert coplanar edges for endpoints of triangles
      for (int i=0;i<3;++i){
        int nindex=triangle_vertices[i]->info();
        if ( nindex >=0 )
          insert_constrained_edges_coplanar_case(nindex,triangulation,id_to_CDT_vh);
      }

      //XSL_TAG_CPL_VERT
      //collect edges incident to a point that is the intersection of two coplanar faces.
      //This ensure that triangulations are compatible.
      if (it_fb!=faces_boundary.end()) //is there at least one intersection point on the boundary of the face?
      {
        for (typename CDT::Finite_vertices_iterator vit=triangulation.finite_vertices_begin(),
                                                vit_end=triangulation.finite_vertices_end();vit_end!=vit;++vit)
        {
          //skip original vertices (that are not nodes) and non-coplanar facet issued vertices
          //(this is working because intersection points between coplanar facets are the first inserted)
          if ( vit->info() < 0 || vit->info() >= number_coplanar_vertices) continue;
          std::map< int,std::set<int> >::iterator res=coplanar_constraints.insert(std::make_pair(vit->info(),std::set<int>())).first;
          //turn around the vertex and get incident edge
          typename CDT::Edge_circulator  start=triangulation.incident_edges(vit);
          typename CDT::Edge_circulator  curr=start;
          do{
            if (triangulation.is_infinite(*curr) ) continue;
            typename CDT::Edge mirror_edge=triangulation.mirror_edge(*curr);
            if ( triangulation.is_infinite( curr->first->vertex(curr->second) ) ||
                 triangulation.is_infinite( mirror_edge.first->vertex(mirror_edge.second) ) )
              continue; //skip edges that are on the boundary of the triangle (these are already constrained)
            //insert edges in the set of constraints
            int nindex =
              curr->first->vertex( (curr->second+1)%3 )==static_cast<typename CDT::Vertex_handle>(vit)?
              (curr->second+2)%3:(curr->second+1)%3;
            typename CDT::Vertex_handle vn=curr->first->vertex(nindex);
            if ( vit->info() > vn->info() ) continue; //take only one out of the two edges + skip negative vn->info()
            CGAL_assertion(vn->info()>=0);
            res->second.insert( vn->info() );
          }while(start!=++curr);
        }

// this is a working alternative that should be slower
//        for (typename CDT::Finite_edges_iterator eit=triangulation.finite_edges_begin(),
//                                             eit_end=triangulation.finite_edges_end();eit_end!=eit;++eit)
//        {
//          typename CDT::Edge mirror_edge=triangulation.mirror_edge(*eit);
//          if ( triangulation.is_infinite( eit->first->vertex(eit->second) ) ||
//              triangulation.is_infinite( mirror_edge.first->vertex(mirror_edge.second) ) )
//            continue; //skip edges that are on the boundary of the triangle (these are already constrained)
//          typename CDT::Vertex_handle v1=eit->first->vertex( (eit->second+1)%3 ),
//                                      v2=eit->first->vertex( (eit->second+2)%3 );
//          if (v1->info()<0 || v2->info()<0) continue;
//          if ( v1->info() > v2->info() ) std::swap(v1,v2);
//          coplanar_constraints.insert(std::make_pair(v1->info(),std::set<int>())).first->second.insert(v2->info());
//        }
      }


      //create a modifier to insert nodes and copy the triangulation of the face
      //inside the polyhedron
      internal_IOP::Triangulate_a_face<typename Polyhedron::HalfedgeDS, NestedFacetConstruct, NewNodeVertexVisitor, PolyhedronPointPMap> modifier(
        f, nodes, node_ids, node_to_polyhedron_vertex, edge_to_hedge, triangulation, facet_construct, node_vertex_visitor, ppmap);

      CGAL_assertion(P->is_valid());
      P->delegate(modifier);
      CGAL_assertion(P->is_valid());

      //3) mark halfedges that are common to two polyhedral surfaces
      //recover halfedges inserted that are on the intersection
      for (std::list<std::pair<int,int> >::iterator it_cst=constrained_edges.begin();it_cst!=constrained_edges.end();++it_cst)
      {
        typename std::map<std::pair<int,int>,Halfedge_handle>::iterator it_poly_hedge=edge_to_hedge.find(*it_cst);
        //we cannot have an assertion here in the case an edge or part of an edge is a constraints.
        //Indeed, the graph_of_constraints report an edge 0,1 and 1,0 for example while only one of the two
        //is defined as one of them defines an adjacent face
        //CGAL_assertion(it_poly_hedge!=edge_to_hedge.end());
        if( it_poly_hedge!=edge_to_hedge.end() ){
          if ( border_halfedges.insert( std::make_pair(Halfedge_const_handle(it_poly_hedge->second),*it_cst) ).second)
          {
            put(m_edge_mark_pmap,std::make_pair(it_poly_hedge->second,P),true);
            put(m_edge_mark_pmap,std::make_pair(it_poly_hedge->second->opposite(),P),true); //setting the opposite is only needed for border edges (done in adjacent triangle otherwise)
          }
          update_edge_per_polyline(P,it_poly_hedge->first,it_poly_hedge->second);
        }
        else{
          //WARNING: in few case this is needed if the marked edge is on the border
          //to optimize it might be better to only use sorted pair. TAG_SLXX1
          std::pair<int,int> opposite_pair(it_cst->second,it_cst->first);
          it_poly_hedge=edge_to_hedge.find(opposite_pair);
          CGAL_assertion( it_poly_hedge!=edge_to_hedge.end() );

          if ( border_halfedges.insert( std::make_pair(Halfedge_const_handle(it_poly_hedge->second),opposite_pair) ).second )
          {
            put(m_edge_mark_pmap,std::make_pair(it_poly_hedge->second,P),true);
            put(m_edge_mark_pmap,std::make_pair(it_poly_hedge->second->opposite(),P),true); //setting the opposite is only needed for border edges (done in adjacent triangle otherwise)
          }
          update_edge_per_polyline(P,it_poly_hedge->first,it_poly_hedge->second);
        }
      }
    }
    output_builder(border_halfedges, nodes, an_edge_per_polyline, is_node_of_degree_one, polyhedron_to_map_node_to_polyhedron_vertex);
  }

  template <class PolylineOfHalfedgeOutputIterator, class Marked_set>
  PolylineOfHalfedgeOutputIterator
  explicitly_compute_polylines(
    Polyhedron* P,
    const Marked_set& is_marked,
    PolylineOfHalfedgeOutputIterator out)
  {
    typedef std::pair< const std::pair<int,int>,
                       std::pair< std::map<Polyhedron*,Halfedge_handle>,
                                  std::pair<bool,int> > >  Complicated_pair;
    BOOST_FOREACH(
      Complicated_pair& p,
      an_edge_per_polyline)
    {
      const std::pair<bool,int>& reversed_and_nbpts = p.second.second;
      Halfedge_handle hedge = p.second.first[P];
      std::vector<Halfedge_handle> polyline;
      int nbsegments=reversed_and_nbpts.second-1;
      polyline.reserve( nbsegments );
      polyline.push_back( reversed_and_nbpts.first?hedge->opposite():hedge );
      for (int i=1; i<nbsegments; ++i)
      {
        hedge = internal_IOP::
          next_marked_halfedge_around_target_vertex (hedge, is_marked);
        polyline.push_back( hedge );
      }
      *out++=polyline;
    }

    return out;
  }

};

}//namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H
