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
// $UR$
// $Id$
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H
#define CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/internal/corefinement/Combinatorial_map_for_corefinement.h> 

#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <boost/optional.hpp>

#include <fstream>
#include <sstream>

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
//  --filtered_order_around_edge: can be done using the original supporting planes
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
namespace CGAL
{
  
  namespace internal_IOP
  {
    template <class Polyhedron>
    struct Compare_unik_address{
      typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
      typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
      typedef typename Polyhedron::Halfedge               Halfedge;
      
      bool operator()(Halfedge_handle h1,Halfedge_handle h2) const {
        Halfedge* ph1=&(*h1) < &(*h1->opposite()) ? &(*h1) : &(*h1->opposite());
        Halfedge* ph2=&(*h2) < &(*h2->opposite()) ? &(*h2) : &(*h2->opposite());
        return  ph1 < ph2; 
      }

      bool operator()(Halfedge_const_handle h1,Halfedge_const_handle h2) const {
        const Halfedge* ph1=&(*h1) < &(*h1->opposite()) ? &(*h1) : &(*h1->opposite());
        const Halfedge* ph2=&(*h2) < &(*h2->opposite()) ? &(*h2) : &(*h2->opposite());
        return  ph1 < ph2; 
      }
    };

  template <class Polyhedron>
  struct Compare_address{
    typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
    typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
    typedef typename Polyhedron::Halfedge               Halfedge;
    
    bool operator()(Halfedge_handle h1,Halfedge_handle h2) const {
      return  &(*h1) < &(*h2); 
    }

    bool operator()(Halfedge_const_handle h1,Halfedge_const_handle h2) const {
      return  &(*h1) < &(*h2); 
    }
  };

  template <class Polyhedron>
  class Non_intersection_halfedge{
    typedef std::map< typename Polyhedron::Halfedge_const_handle,
                      std::pair<int,int>,
                      Compare_unik_address<Polyhedron> 
                    >  Intersection_hedges_set;
    Intersection_hedges_set intersection_hedges_;
  public:  
    Non_intersection_halfedge(const Intersection_hedges_set& the_set) : intersection_hedges_(the_set){}
  
  
    bool operator()(typename Polyhedron::Halfedge_const_handle h) const
    {
      return intersection_hedges_.find(h)==intersection_hedges_.end();
    }
  };

    
  template <class HDS>
  class Triangulate_a_face : public CGAL::Modifier_base<HDS> {
    typedef typename HDS::Halfedge_handle Halfedge_handle;
    typedef typename HDS::Vertex_handle   Vertex_handle;
    typedef typename HDS::Face_handle     Face_handle;
    typedef typename HDS::Vertex          Vertex;
    typedef typename HDS::Halfedge        Halfedge;
    typedef typename HDS::Face            Face;
    
    //data members
    Face_handle current_face;
    std::map<int,typename Vertex::Point >                  nodes_;
    std::map<int,Vertex_handle>&                           node_to_polyhedron_vertex_;
    std::map<std::pair<int,int>,Halfedge_handle>&          edge_to_hedge_;
    std::vector<std::pair<int,int> >                       edges_to_create_;
    std::vector<CGAL::cpp0x::tuple<int,int,int> >          faces_to_create_;
    
    typename HDS::Halfedge::Base*
    unlock_halfedge(Halfedge_handle h){
      return static_cast<typename HDS::Halfedge::Base*>(&(*h));
    }
    
    typename HDS::Face::Base*
    unlock_face(Face_handle f){
      return static_cast<typename HDS::Face::Base*>(&(*f));
    }    
    
  public:
    
    template <class Node_vector,class Triangulation>
    Triangulate_a_face( Face_handle face,
                        const Node_vector& nodes,
                        const std::vector<int>& node_ids,
                        std::map<int,Vertex_handle>& node_to_polyhedron_vertex,
                        std::map<std::pair<int,int>,Halfedge_handle>& edge_to_hedge,
                        const Triangulation& triangulation)
    :current_face(face),node_to_polyhedron_vertex_(node_to_polyhedron_vertex),edge_to_hedge_(edge_to_hedge)
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
        faces_to_create_.push_back( CGAL::cpp0x::make_tuple( v0->info(),v1->info(),v2->info() ) );
      }
    }


  
    void operator()( HDS& hds) {
//      std::cerr << "node_to_polyhedron_vertex_"<< std::endl;
//      for (typename std::map<int,Vertex_handle>::iterator it=node_to_polyhedron_vertex_.begin();it!=node_to_polyhedron_vertex_.end();++it)
//        std::cerr << it->first << " " << &(*(it->second)) << std::endl;
      
      //insert the intersection point interior to the face inside the polyhedron and
      //save their Polyhedron::vertex_handle 
      for (typename std::map<int,typename Vertex::Point>::iterator it=nodes_.begin();it!=nodes_.end();++it)
      {
        Vertex_handle v=hds.vertices_push_back(Vertex(it->second));
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
      
      std::vector<CGAL::cpp0x::tuple<int,int,int> >::iterator it=faces_to_create_.begin();
      
      //create the new faces and update adjacencies
      while (true)
      {
        int i=cpp0x::get<0>(*it),j=cpp0x::get<1>(*it),k=cpp0x::get<2>(*it);
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
          current_face=hds.faces_push_back(Face());
        else
          break;
      }
    }
  };    
 
  } //namespace internal_IOP

  
//Considering the plane with normal vector [O_prime,O] and containing O. 
//We define the counterclockwise order around O when looking from the side of the plane 
//into which the vector [O_prime,O] is pointing.
//We consider the portion of the plane defined by rotating a ray starting at O
//from the planar projection of P1 to the planar projection of P2 in counterclockwise order.
//The predicates indicates whether the planar projection of point Q lies in this portion of the plane.
//Preconditions:
//  O_prime,O,P1 are not collinear
//  O_prime,O,P2 are not collinear
//  O_prime,O,Q are not collinear
//  O_prime,O,P1,Q are not coplanar or coplanar_orientation(O,O_prime,P1,Q)==NEGATIVE
//  O_prime,O,P2,Q are not coplanar or coplanar_orientation(O,O_prime,P2,Q)==NEGATIVE
template <class Kernel>
bool  is_in_interior_of_object(
    const typename Kernel::Point_3& O_prime,const typename Kernel::Point_3& O,
    const typename Kernel::Point_3& P1,const typename Kernel::Point_3& P2,
    const typename Kernel::Point_3& Q)
{
  //guarantee to have non-flat triangles
  CGAL_precondition( !collinear(O_prime,O,P1) );
  CGAL_precondition( !collinear(O_prime,O,P2) );
  CGAL_precondition( !collinear(O_prime,O,Q)  );

  //no two triangles are coplanar and on the same side of their common edge
  CGAL_precondition( !coplanar(O_prime,O,P1,Q) || coplanar_orientation(O,O_prime,P1,Q)==NEGATIVE );
  CGAL_precondition( !coplanar(O_prime,O,P2,Q) || coplanar_orientation(O,O_prime,P2,Q)==NEGATIVE );
  
  Sign s0 = sign( determinant(O-O_prime,P1-O,P2-O) );
  
  if ( s0==ZERO ){
    //O, O_prime, P1 and P2 are coplanar
    Orientation o=orientation(O_prime,O,P1,Q);
    CGAL_precondition(o!=COPLANAR);
    return o==POSITIVE;
  }
  
  //O, O_prime, P1 and P2 are not coplanar
  Sign s1 = sign( determinant(O-O_prime,P1-O,Q -O) );
  Sign s2 = sign( determinant(O-O_prime,Q -O,P2-O) );

  if (s0 == POSITIVE) // the angle P1,O,P2 is smaller that Pi.
    return ( s1 == POSITIVE ) && ( s2 ==POSITIVE ); //true if the angles P1,O,Q and Q,O,P2 are smaller than Pi
  else
    return ( s1 != NEGATIVE ) || ( s2 != NEGATIVE ); //true if the angle P1,O,Q or the angle Q,O,P2 is smaller than or equal to Pi
}

//import into the combinatorial map facets in the given range.
//they are supposed to be in the same connected component.
//two volume are created (each facets gives two opposite orientation 2-cell in the map)
template<class Polyhedron, class Map, class Face_iterator, class Non_special_edge_predicate,class Halfedge_to_dart_map_ >
typename Map::Dart_handle import_from_polyhedron_subset(  Map& amap,
                                                          Face_iterator faces_begin,
                                                          Face_iterator faces_end,
                                                          const Non_special_edge_predicate& is_non_special_edge,
                                                          Halfedge_to_dart_map_& selected_hedge_to_dart,
                                                          int mark_index
  )
{
  typedef typename Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
  typedef std::map < Halfedge_const_handle, typename Map::Dart_handle,internal_IOP::Compare_address<Polyhedron> > Halfedge_to_dart_map;
   
  Halfedge_to_dart_map hedge_to_dart;
  typename Map::Dart_handle first_dart = NULL;
  // First traversal to build the darts and link them.
  for (Face_iterator it_face = faces_begin; it_face != faces_end; ++it_face)
  {
    Halfedge_const_handle start=(*it_face)->halfedge();
    
    CGAL_precondition(start->next()!=start);
    
    Halfedge_const_handle current=start;
    typename Map::Dart_handle prev = NULL;
    typename Map::Dart_handle first_dart_of_face = NULL;
    do
    {
      typename Map::Dart_handle d = amap.create_dart();
      amap.template link_beta<3>(d,amap.create_dart()); //for opposite volume
      hedge_to_dart[current] = d;
             
      if (prev != NULL){
        amap.template link_beta<1>(prev, d);
        amap.template link_beta<1>(d->beta(3),prev->beta(3));//for opposite volume
      }
      else 
      {
        first_dart_of_face = d;
        if (first_dart==NULL) first_dart=d;
      }
      
      if ( is_non_special_edge (current) ){
        if ( !current->is_border_edge() ){
          CGAL_assertion(current != current->opposite());
          typename Halfedge_to_dart_map::iterator it = hedge_to_dart.find(current->opposite());
          if (it != hedge_to_dart.end()){ //link the opposites halfedges only when both corresponding darts have been created
            amap.template link_beta<2>(d, it->second);
            amap.template link_beta<2>(d->beta(3), it->second->beta(3));//for opposite volume
          }
        }
      }
      else{
        typename Halfedge_to_dart_map_::iterator it_hedge_map=selected_hedge_to_dart.find(current);
         //all marked hedges are not the selected one for its polyline
        if ( it_hedge_map!=selected_hedge_to_dart.end() ) it_hedge_map->second=d;
        //darts d and d->beta(3) are special edges
        amap.mark(d,mark_index);
        amap.mark(d->beta(3),mark_index);
      }
      prev = d;
      current=current->next();
    }
    while (current != start);
    amap.template link_beta<1>(prev, first_dart_of_face);
    amap.template link_beta<1>(first_dart_of_face->beta(3),prev->beta(3));//for opposite volume
  }

  // Second traversal to update the geometry.
  // We run one again through the facets of the HDS.
  for (Face_iterator it_face = faces_begin; it_face != faces_end; ++it_face)
  {
    Halfedge_const_handle start=(*it_face)->halfedge();
    Halfedge_const_handle current=start;
    do
    {
      typename Map::Dart_handle d = hedge_to_dart[current]; // Get the dart associated to the Halfedge
      if (d->template attribute<0>() == NULL)
      {	    
        amap.template set_attribute<0>(d,
           amap.template create_attribute<0>(current->opposite()->vertex()->point()));
      }
      current=current->next();
    }
    while (current != start);
  }
  
  return first_dart;
}

 //turn around the target vertex of dart to find a marked dart
template <class Combinatorial_map_3>
boost::optional<typename Combinatorial_map_3::Dart_handle> 
next_marked_dart_around_target_vertex(
  const Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  int mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=dart->beta(1);
  while ( ! final_map.is_marked(next,mark_index) ){
    if (next->is_free(2) )//we reach a boundary
      return  boost::optional<typename Combinatorial_map_3::Dart_handle>();
    next=next->beta(2)->beta(1);
  }
  if (next == dart) //no new dart have been found  
    return  boost::optional<typename Combinatorial_map_3::Dart_handle>();
  CGAL_precondition(&dart->beta(1)->template attribute<0>()->point() == &next->template attribute<0>()->point());
  return boost::optional<typename Combinatorial_map_3::Dart_handle> (next);
}

//turn around the target vertex of dart to find a marked dart
//with expected_target as target vertex
template <class Combinatorial_map_3>
typename Combinatorial_map_3::Dart_handle
get_next_marked_dart_around_target_vertex(
  const Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  int mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=dart->beta(1);
  while ( !final_map.is_marked(next,mark_index) ){
    CGAL_assertion( !next->is_free(2) );
    next=next->beta(2)->beta(1);
    CGAL_assertion(next != dart);
  }
  CGAL_precondition(&dart->beta(1)->template attribute<0>()->point() == &next->template attribute<0>()->point());
  return next;
}

//turn around the source vertex of dart to find a marked dart
//with expected_source as source vertex
template <class Combinatorial_map_3>
typename Combinatorial_map_3::Dart_handle 
get_next_marked_dart_around_source_vertex(  
  const Combinatorial_map_3& final_map,
  typename Combinatorial_map_3::Dart_handle dart,
  int mark_index)
{
  CGAL_precondition(final_map.is_marked(dart,mark_index));
  typename Combinatorial_map_3::Dart_handle next=dart->beta(0);
  while ( ! final_map.is_marked(next,mark_index) ){ 
    CGAL_assertion( !next->is_free(2) );
    next=next->beta(2)->beta(0);
    CGAL_assertion(next != dart);
  }
  CGAL_precondition(&dart->template attribute<0>()->point() == &next->beta(1)->template attribute<0>()->point());
  return next;
}

//given two marked darts, this function links these two darts with beta<2>
//but in addition it follows the marked darts connected to the same vertex
//(there should be only one) to connect them all together
//( this function is a kind of zipper ;) )
template <class Combinatorial_map_3,class Node_vector>
void sew_2_marked_darts( Combinatorial_map_3& final_map,
                         typename Combinatorial_map_3::Dart_handle dart_1 , 
                         typename Combinatorial_map_3::Dart_handle dart_2 ,
                         int mark_index,
                         const Node_vector& nodes,
                         const std::pair<int,int>& indices,
                         const std::pair<bool,int>& polyline_info)
{
  CGAL_precondition( dart_1->is_free(2) );
  CGAL_precondition( dart_2->is_free(2) );
  CGAL_precondition( final_map.is_marked(dart_1,mark_index) );
  CGAL_precondition( final_map.is_marked(dart_2,mark_index) );
  CGAL_precondition( dart_1->template attribute<0>()->point() == dart_2->beta(1)->template attribute<0>()->point() );
  CGAL_precondition( dart_1->beta(1)->template attribute<0>()->point() == dart_2->template attribute<0>()->point() );
  
  int src_index = ( ( indices.first < indices.second) ==  polyline_info.first )
                  ? indices.second:indices.first;
  
  if ( dart_1->template attribute<0>()->point() != nodes[ src_index ] ) std::swap(dart_1,dart_2);
  
  int nb_segs=polyline_info.second-1,k=1;
  
  do{
    CGAL_precondition( final_map.template is_sewable<2>(dart_1,dart_2) );
    final_map.template sew<2>(dart_1,dart_2);
    
    if (k==nb_segs) break;
      
    dart_1=get_next_marked_dart_around_target_vertex(final_map,dart_1,mark_index);
    dart_2=get_next_marked_dart_around_source_vertex(final_map,dart_2,mark_index);
  }
  while(++k);
}

//not_top and not_down are two darts from volumes that get merged with an existing
//other one because of a set of identical coplanar triangles.
//top and down is the dart of the volumes "replacing" that of not_top and not down respectively, 
//The function is considering all triangles that are bounded by a cycle of marked edges.
//The volume not_top and not_down are part of are those that will disappear at the
//end of the main algorithm.
//( this function is a kind of facet gluer ;) )
template <class Combinatorial_map_3>
void sew_3_marked_darts( Combinatorial_map_3& final_map,
                         typename Combinatorial_map_3::Dart_handle not_top , 
                         typename Combinatorial_map_3::Dart_handle not_down ,
                         typename Combinatorial_map_3::Dart_handle top , 
                         typename Combinatorial_map_3::Dart_handle down ,
                         int mark_index,
                         std::set<typename Combinatorial_map_3::Dart_handle>& darts_to_remove)
{
  typedef boost::optional<typename Combinatorial_map_3::Dart_handle> O_Dart_handle;

  if ( not_top->template attribute<3>()->info().is_empty ){
    CGAL_assertion(not_down->template attribute<3>()->info().is_empty);
    return;
  }
   
  CGAL_assertion(!not_down->template attribute<3>()->info().is_empty);

  //merge attribute of the two volumes:
  internal_IOP::Volume_on_merge merge_attributes;
  merge_attributes(*top->template attribute<3>(),*not_top->template attribute<3>());
  merge_attributes(*down->template attribute<3>(),*not_down->template attribute<3>());
  
  //set volume attributes as empty to avoid double sew_3 of the same topological disk of triangles
  not_top->template attribute<3>()->info().is_empty=true;
  not_down->template attribute<3>()->info().is_empty=true;
  
  CGAL_precondition( final_map.is_marked(not_top,mark_index) && final_map.is_marked(top,mark_index) );
  CGAL_precondition( final_map.is_marked(not_down,mark_index) && final_map.is_marked(down,mark_index) );
  CGAL_precondition( not_top->template attribute<0>()->point() == not_down->beta(1)->template attribute<0>()->point() );
  CGAL_precondition( not_top->beta(1)->template attribute<0>()->point() == not_down->template attribute<0>()->point() );
  CGAL_precondition( not_top->template attribute<0>()->point() == top->template attribute<0>()->point() );
  CGAL_precondition( not_down->template attribute<0>()->point() == down->template attribute<0>()->point() );
  
  CGAL_assertion( top->beta(3)==down );

  //set to be removed the darts of the two no longer used volumes
  typename Combinatorial_map_3::Dart_handle start=not_top;
  do
  {
    CGAL_assertion(!not_top->is_free(3));
    darts_to_remove.insert(not_top);   darts_to_remove.insert(not_top->beta(1)); darts_to_remove.insert(not_top->beta(1)->beta(1));
    darts_to_remove.insert(not_top->beta(3));   darts_to_remove.insert(not_top->beta(3)->beta(1)); darts_to_remove.insert(not_top->beta(3)->beta(1)->beta(1));
    O_Dart_handle current_1=next_marked_dart_around_target_vertex(final_map,not_top,mark_index);
    CGAL_precondition(current_1);
    not_top=*current_1;
  }
  while(not_top!=start);
}

template<class Tag>
struct Halfedge_marker{
  template <class Halfedge_handle>
  static void mark(Halfedge_handle){}
};

template<>
struct Halfedge_marker<Tag_true>{
  template <class Halfedge_handle>
  static void mark(Halfedge_handle h){h->set_mark();}
};

template<class Polyhedron,class Kernel=typename Polyhedron::Traits::Kernel,class Mark_intersection_halfedges=Tag_false>
class Node_visitor_refine_polyhedra{
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
  typedef CGAL::Exact_predicates_exact_constructions_kernel             Exact_kernel;
  typedef CGAL::Triangulation_vertex_base_with_info_2<int,Exact_kernel> Vbi;
  typedef CGAL::Constrained_triangulation_face_base_2<Exact_kernel>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi,Fb>                  TDS_2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Exact_kernel,TDS_2,CGAL::No_intersection_tag> CDT;  //DO WE NEED DELAUNAY????  
  #endif
  
  typedef std::map<Halfedge_handle,Polyhedron*,Cmp_unik_ad>            Hedge_to_polyhedron_map;

  typedef std::vector<int>                                             Node_ids;
  typedef std::set<int>                                                Node_id_set; //avoid having duplicated node on edge of coplanar triangles
  typedef std::map< Face_handle,Node_ids,Cmp_handle >                  In_face_map; 
  typedef std::map< Halfedge_handle,Node_id_set,Cmp_unik_ad >          In_halfedge_map;
  //to keep the correspondance between node_id and vertex_handle in each polyhedron
  typedef std::map<int,Vertex_handle> Node_to_polyhedron_vertex_map;
  typedef std::map<Polyhedron*, Node_to_polyhedron_vertex_map > Poly_to_map_node;
  //to maintain an polyhedron halfedge on each polyline + pair<bool,int>
  //with first = "is the key (pair<int,int>) was reversed?" and second is the number of edges +1 in the polyline
  typedef std::map< std::pair<int,int>, std::pair< std::map<Polyhedron*,Halfedge_handle>,std::pair<bool,int> > > An_edge_per_polyline_map;  
  //to handle coplanar halfedge of polyhedra that are full in the intersection
  typedef std::map< int,Halfedge_handle >                              Node_to_target_of_hedge_map;
  typedef std::map< Polyhedron*,Node_to_target_of_hedge_map>           Poly_to_vertices_on_intersection_map;

  //Combinatorial map typedefs
  typedef internal_IOP::Item_with_points_and_volume_info<Kernel,Polyhedron> Items;
  typedef CGAL::Combinatorial_map<3,Items> Combinatorial_map_3_;
  typedef typename Combinatorial_map_3_::Dart_handle Dart_handle;
  
//data members
  Hedge_to_polyhedron_map               hedge_to_polyhedron;
  In_face_map                           in_face;
  In_halfedge_map                       in_hedge;
  std::map< int,std::set<int> >         graph_of_constraints;
  std::map< int,std::set<int> >         coplanar_constraints;
  An_edge_per_polyline_map              an_edge_per_polyline;
  typename An_edge_per_polyline_map::iterator last_polyline;
  Poly_to_vertices_on_intersection_map  poly_to_vertices_on_inter;
  Poly_to_map_node                      polyhedron_to_map_node_to_polyhedron_vertex;
  std::set<int>                         non_manifold_nodes; //contain nodes that are original vertices of input polyhedron and that neighborhood is not a topological disk
  std::map<Vertex_handle,int>           nodes_that_are_original_vertices;//to keep the correspondance between original polyhedron vertices that are also nodes
  
  Combinatorial_map_3_*                 final_map_ptr;
  Combinatorial_map_3_&                 final_map() {return *final_map_ptr;}
  bool final_map_comes_from_outside;
  //   new_hedge    hedge
  //  ----------->   ----------->
  //               v
  //  <-----------   <-----------
  //   new_opposite     opposite 
  //  
  Vertex_handle split_edge( Halfedge_handle hedge,
                            const typename Kernel::Point_3& point,
                            Polyhedron& P)
  {
    internal_IOP::Split_halfedge_at_point<typename Polyhedron::HalfedgeDS> delegated(hedge,point);
    P.delegate( delegated );
    CGAL_assertion(P.is_valid());
    
    Vertex_handle v=boost::prior(P.vertices_end());
    CGAL_assertion(v->point()==point);
    return v;
  }

  //sort node ids so that we can split the hedge
  //consecutively
  template <class Node_vector>
  void sort_vertices_along_hedge(std::vector<int>& node_ids,Halfedge_handle hedge,const Node_vector& nodes)
  {
    std::sort(node_ids.begin(),
              node_ids.end(),
              internal_IOP::Order_along_a_halfedge<Polyhedron,Node_vector,Is_polyhedron_const>(hedge,nodes)
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
      #ifndef NDEBUG
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

  int node_index_of_incident_vertex(Halfedge_const_handle h,
                                    const std::map<Halfedge_const_handle,std::pair<int,int>,Cmp_unik_ad >& border_halfedges)
  {
    //WARNING this may be expensive
    typedef std::map<Halfedge_const_handle,std::pair<int,int>,Cmp_unik_ad > Border_halfedges_map;

    Halfedge_const_handle start=h;
    Halfedge_const_handle curr=start;
    do{
      typename Border_halfedges_map::const_iterator it_border=border_halfedges.find( curr );
      if (it_border!=border_halfedges.end())
        return it_border->first==curr?it_border->second.second:it_border->second.first;      
      curr=curr->next()->opposite();
    }while(curr!=start);
    
    return -1;
  }
  
  template <class Node_vector>
  bool filtered_order_around_edge(int O_prime_index,
                                  int O_index,
                                  int P1_index,
                                  int P2_index,
                                  int Q_index,
                                  Vertex_handle P1,
                                  Vertex_handle P2,
                                  Vertex_handle Q,
                                  const Node_vector& nodes)
  {
    try{
      return is_in_interior_of_object<typename Node_vector::Ikernel>(
        nodes.interval_node(O_prime_index),
        nodes.interval_node(O_index),
        P1_index == -1 ? nodes.to_interval(P1->point()): nodes.interval_node(P1_index),
        P2_index == -1 ? nodes.to_interval(P2->point()): nodes.interval_node(P2_index),
        Q_index  == -1 ? nodes.to_interval(Q->point()) : nodes.interval_node(Q_index )
      );
    }
    catch(Uncertain_conversion_exception&){
      return is_in_interior_of_object<typename Node_vector::Exact_kernel>(
        nodes.exact_node(O_prime_index),
        nodes.exact_node(O_index),
        P1_index == -1 ? nodes.to_exact(P1->point()): nodes.exact_node(P1_index),
        P2_index == -1 ? nodes.to_exact(P2->point()): nodes.exact_node(P2_index),
        Q_index  == -1 ? nodes.to_exact(Q->point()) : nodes.exact_node(Q_index )
      );
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

//======================================================================//
//functions internally used to glue piece of the final combinatorial map//
//======================================================================//


  //-----first polyhedron
template <class Halfedge_to_dart_map>
inline Dart_handle get_associated_dart(Halfedge_handle hedge,Halfedge_to_dart_map& selected_hedge_to_dart){
  typename Halfedge_to_dart_map::iterator it_saved_dart=selected_hedge_to_dart.find(hedge);
  CGAL_assertion(it_saved_dart!=selected_hedge_to_dart.end());
  return it_saved_dart->second;
}

//first_hedge defines four volumes, second_hedge only two
//first_poly is not needed as inside/outside volume is update during the merge
//of the sew. Only second_poly is needed
template <class Node_vector,class Border_halfedges_map,class Halfedge_to_dart_map>
void sew_2_three_volumes_case(  Halfedge_handle first_hedge, Halfedge_handle second_hedge,
                                const std::pair<int,int>& indices,
                                const Node_vector& nodes,
                                Border_halfedges_map& border_halfedges,
                                Halfedge_to_dart_map& selected_hedge_to_dart, 
                                Polyhedron* /*first_poly*/, Polyhedron* second_poly,
                                int mark_index,
                                std::set<Dart_handle>& darts_to_remove,
                                const std::pair<bool,int>& polyline_info)
{
  bool took_opposite=second_hedge->is_border();
  if (took_opposite) second_hedge=second_hedge->opposite();
  
  Vertex_handle P1=first_hedge->opposite()->next()->vertex();
  Vertex_handle P2=first_hedge->next()->vertex();
  //    when looking from the side of indices.second, the interior of the first polyhedron is described 
  //    by turning counterclockwise from P1 to P2
  
  Vertex_handle Q = second_hedge->next()->vertex();
  
  //check if the third point of each triangular face is an original point (stay -1)
  //or a intersection point (in that case we need the index of the corresponding node to
  //have the exact value of the point)      
  int index_p1=node_index_of_incident_vertex(first_hedge->opposite()->next(),border_halfedges);
  int index_p2=node_index_of_incident_vertex(first_hedge->next(),border_halfedges);
  int index_q =node_index_of_incident_vertex(second_hedge->next(),border_halfedges);      

  //Recover the dart that will be the start point of the different sewing
  //  dof_X_outside = dart of face of , meaning the triangle containing the
  //  point X and part of the volume outside of the corresponding polyhedron      
  //-----first polyhedron
  Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),selected_hedge_to_dart);
  Dart_handle dof_P2_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
  //-----second polyhedron
  Dart_handle dof_Q_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);
  
  if (index_p1!=-1 && index_p1==index_q){
    Dart_handle top=dof_P1_outside->beta(3), not_top=took_opposite?dof_Q_outside->beta(3):dof_Q_outside;
    Dart_handle down=dof_P1_outside, not_down=took_opposite?dof_Q_outside:dof_Q_outside->beta(3);

    if ( top->template attribute<3>()->info().is_empty ) std::swap(not_top,top);
    if ( down->template attribute<3>()->info().is_empty ) std::swap(not_down,down);
    CGAL_assertion( !top->template attribute<3>()->info().is_empty );
    CGAL_assertion( !down->template attribute<3>()->info().is_empty );

    sew_2_marked_darts( final_map(),top             , dof_P2_outside->beta(3)     ,mark_index, nodes, indices, polyline_info); //P1P2 or QP2      
    sew_2_marked_darts( final_map(),dof_P2_outside  , down                        ,mark_index, nodes, indices, polyline_info); //P2Q or P2P1
    sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,darts_to_remove);
    
    return;
  }

  if (index_p2!=-1 && index_p2==index_q){
    Dart_handle top=dof_P2_outside->beta(3), not_top=took_opposite?dof_Q_outside:dof_Q_outside->beta(3);
    Dart_handle down=dof_P2_outside, not_down=took_opposite?dof_Q_outside->beta(3):dof_Q_outside;

    if ( top->template attribute<3>()->info().is_empty ) std::swap(not_top,top);
    if ( down->template attribute<3>()->info().is_empty ) std::swap(not_down,down);
    CGAL_assertion( !top->template attribute<3>()->info().is_empty );
    CGAL_assertion( !down->template attribute<3>()->info().is_empty );

    sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , top              ,mark_index, nodes, indices, polyline_info); //P1Q or P1P2      
    sew_2_marked_darts( final_map(),down                    , dof_P1_outside   ,mark_index, nodes, indices, polyline_info); //QP1 or P2P1
    sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,darts_to_remove);
    
    return;
  }
  
  bool Q_is_between_P1P2 = filtered_order_around_edge(indices.first,indices.second,index_p1,index_p2,index_q,P1,P2,Q,nodes);


  if (Q_is_between_P1P2)
  {
    // poly_first  - poly_second            = took_opposite?P1Q:QP2
    // poly_second - poly_first             = {0}
    // poly_first \cap poly_second          = took_opposite?QP2:P1Q
    // opposite( poly_first U poly_second ) = P2P1
    sew_2_marked_darts( final_map(),dof_P1_outside->beta(3)                            , took_opposite?dof_Q_outside:dof_Q_outside->beta(3)  ,mark_index, nodes, indices, polyline_info); //P1Q
    sew_2_marked_darts( final_map(),took_opposite?dof_Q_outside->beta(3):dof_Q_outside , dof_P2_outside->beta(3)                             ,mark_index, nodes, indices, polyline_info); //QP2
    sew_2_marked_darts( final_map(),dof_P2_outside                                     , dof_P1_outside                                      ,mark_index, nodes, indices, polyline_info); //P2P1
    dof_P1_outside->template attribute<3>()->info().outside.insert(second_poly); //update P2P1 outside poly
  }
  else
  {
    // poly_first  - poly_second            = P1P2
    // poly_second - poly_first             = took_opposite?QP1:P2Q
    // poly_first \cap poly_second          = {0}
    // opposite( poly_first U poly_second ) = took_opposite?P2Q:QP1
    sew_2_marked_darts( final_map(),dof_P2_outside                                     , took_opposite?dof_Q_outside:dof_Q_outside->beta(3)  ,mark_index, nodes, indices, polyline_info); //P2Q
    sew_2_marked_darts( final_map(),took_opposite?dof_Q_outside->beta(3):dof_Q_outside , dof_P1_outside                                      ,mark_index, nodes, indices, polyline_info); //QP1
    sew_2_marked_darts( final_map(),dof_P1_outside->beta(3)                            , dof_P2_outside->beta(3)                             ,mark_index, nodes, indices, polyline_info); //P1P2
    dof_P1_outside->beta(3)->template attribute<3>()->info().outside.insert(second_poly); //update P1P2 outside poly
  }
}

//first_hedge defines two volumes, second_hedge only two
template <class Halfedge_to_dart_map,class Border_halfedges_map,class Node_vector>
void sew_2_two_volumes_case(  Halfedge_handle first_hedge, Halfedge_handle second_hedge,
                              Border_halfedges_map& border_halfedges,
                              Halfedge_to_dart_map& selected_hedge_to_dart,
                              int mark_index,
                              std::set<Dart_handle>& darts_to_remove,
                              const Node_vector& nodes,
                              const std::pair<int,int>& indices,
                              const std::pair<bool,int>& polyline_info)
{
  bool first_took_opposite=first_hedge->is_border();
  if (first_took_opposite) first_hedge=first_hedge->opposite();
  bool second_took_opposite=second_hedge->is_border();
  if (second_took_opposite) second_hedge=second_hedge->opposite();
  
    //-----first polyhedron
  Dart_handle dof_P_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
  //-----second polyhedron
  Dart_handle dof_Q_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);
  
  
  
  
  int index_p =node_index_of_incident_vertex(first_hedge->next(),border_halfedges);
  int index_q =node_index_of_incident_vertex(second_hedge->next(),border_halfedges);

  if (index_p!=-1 && index_q!=-1 && index_p==index_q){
    Dart_handle top=dof_P_outside, not_top=dof_Q_outside->beta(3);
    Dart_handle down=dof_P_outside->beta(3), not_down=dof_Q_outside;
    
    if (first_took_opposite==second_took_opposite)
    {
      top=dof_P_outside->beta(3); not_top=dof_Q_outside->beta(3);
      down=dof_P_outside; not_down=dof_Q_outside;
    }
    
    if ( top->template attribute<3>()->info().is_empty ) std::swap(not_top,top);
    if ( down->template attribute<3>()->info().is_empty ) std::swap(not_down,down);
    CGAL_assertion( !top->template attribute<3>()->info().is_empty );
    CGAL_assertion( !down->template attribute<3>()->info().is_empty );

    sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,darts_to_remove);
    
    return;
  }
  
  
  
  //since the edge is shared, the inside of each polyhedron must be on opposite orientation halfedges
  if (first_took_opposite==second_took_opposite)
  {
    //sew out with in
    sew_2_marked_darts( final_map(),dof_P_outside->beta(3)  , dof_Q_outside  ,mark_index, nodes, indices, polyline_info); //PQ
    sew_2_marked_darts( final_map(),dof_Q_outside->beta(3)  , dof_P_outside  ,mark_index, nodes, indices, polyline_info); //QP
  }
  else
  {
    //sew in with in
  sew_2_marked_darts( final_map(),dof_P_outside         , dof_Q_outside          ,mark_index, nodes, indices, polyline_info); //PQ
  sew_2_marked_darts( final_map(),dof_Q_outside->beta(3), dof_P_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //QP
  }
}

//4 volume case with 2 identical volume
//Q2 is supposed to be identical to P2
template <class Node_vector,class Halfedge_to_dart_map>
void sew_2_four_volumes_case_1(  Halfedge_handle first_hedge, Halfedge_handle second_hedge,
                                 const std::pair<int,int>& indices,
                                 const Node_vector& nodes,
                                 int index_p1, int index_p2, int index_q1,
                                 Halfedge_to_dart_map& selected_hedge_to_dart,
                                 int mark_index,
                                 std::set<Dart_handle>& darts_to_remove,
                                 const std::pair<bool,int>& polyline_info,
                                 bool swap_in_out_Q=false)
{
  Vertex_handle P1=first_hedge->opposite()->next()->vertex();
  Vertex_handle P2=first_hedge->next()->vertex();
  //    when looking from the side of indices.second, the interior of the first polyhedron is described 
  //    by turning counterclockwise from P1 to P2
  Vertex_handle Q1=second_hedge->opposite()->next()->vertex();
  // Vertex_handle Q2=second_hedge->next()->vertex();  
  bool Q1_is_between_P1P2 = filtered_order_around_edge(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes);
  
  
  //Recover the dart that will be the start point of the different sewing
  //  dof_X_outside = dart of face of , meaning the triangle containing the
  //  point X and part of the volume outside of the corresponding polyhedron      
  //-----first polyhedron
  Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),selected_hedge_to_dart);
  Dart_handle dof_P2_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
  //-----second polyhedron
  Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),selected_hedge_to_dart);
  Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);  
  
  if( swap_in_out_Q ){
    dof_Q1_outside=dof_Q1_outside->beta(3);
    dof_Q2_outside=dof_Q2_outside->beta(3);
  }
  
  if (Q1_is_between_P1P2){
      Dart_handle top=dof_Q2_outside->beta(3), not_top=dof_P2_outside->beta(3);
      Dart_handle down=dof_Q2_outside, not_down=dof_P2_outside;
      if ( top->template attribute<3>()->info().is_empty ) std::swap(not_top,top);
      if ( down->template attribute<3>()->info().is_empty ) std::swap(not_down,down);
      CGAL_assertion( !top->template attribute<3>()->info().is_empty );
      CGAL_assertion( !down->template attribute<3>()->info().is_empty );
    
      // poly_first  - poly_second            = P1Q1
      // poly_second - poly_first             = {0}
      // poly_first \cap poly_second          = Q1P2 or Q1Q2
      // opposite( poly_first U poly_second ) = Q2P1 or P2P1
      sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //P1Q1      
      sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , top                     ,mark_index, nodes, indices, polyline_info); //Q1P2 or Q1Q2
      sew_2_marked_darts( final_map(),down                    , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q2P1 or P2P1
      sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,darts_to_remove);

  }
  else{
      Dart_handle top=dof_Q2_outside->beta(3), not_top=dof_P2_outside->beta(3);
      Dart_handle down=dof_Q2_outside, not_down=dof_P2_outside;
      if ( top->template attribute<3>()->info().is_empty ) std::swap(not_top,top);
      if ( down->template attribute<3>()->info().is_empty ) std::swap(not_down,down);
      CGAL_assertion( !top->template attribute<3>()->info().is_empty );
      CGAL_assertion( !down->template attribute<3>()->info().is_empty );

      // poly_first  - poly_second            = {0}
      // poly_second - poly_first             = Q1P1
      // poly_first \cap poly_second          = P1P2 or P1Q2
      // opposite( poly_first U poly_second ) = Q2Q1 or P2Q1
      sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q1P1
      sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , top                     ,mark_index, nodes, indices, polyline_info); //P1P2 or P1Q2
      sew_2_marked_darts( final_map(),down                    , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //Q2Q1 or P2Q1
      sew_3_marked_darts( final_map(),not_top,not_down,top,down,mark_index,darts_to_remove);
  }  
}

template <class Node_vector,class Halfedge_to_dart_map>
bool coplanar_triangles_case_handled(Halfedge_handle first_hedge,Halfedge_handle second_hedge,
                                     const std::pair<int,int>& indices,
                                     const Node_vector& nodes,
                                     int index_p1, int index_p2, int index_q1,int index_q2,
                                     Halfedge_to_dart_map& selected_hedge_to_dart,
                                     int mark_index,
                                     std::set<Dart_handle>& darts_to_remove,
                                     const std::pair<bool,int>& polyline_info)
{
  if( index_p1!=-1 ){
    if (index_p1==index_q1){
      if(index_p2!=-1){
        CGAL_assertion(index_p2!=index_q1);
        if(index_p2==index_q2){
          //-----first polyhedron
          Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_P2_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
          //-----second polyhedron
          Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);      
    
          Dart_handle top_1=dof_P1_outside->beta(3), not_top_1=dof_Q1_outside->beta(3);
          Dart_handle top_2=dof_P2_outside->beta(3), not_top_2=dof_Q2_outside->beta(3);
          Dart_handle down_1=dof_P1_outside, not_down_1=dof_Q1_outside;
          Dart_handle down_2=dof_P2_outside, not_down_2=dof_Q2_outside;
          if ( top_1->template attribute<3>()->info().is_empty ) std::swap(top_1,not_top_1);
          if ( top_2->template attribute<3>()->info().is_empty ) std::swap(top_2,not_top_2);
          if ( down_1->template attribute<3>()->info().is_empty ) std::swap(down_1,not_down_1);
          if ( down_2->template attribute<3>()->info().is_empty ) std::swap(down_2,not_down_2);
          CGAL_assertion( !top_1->template attribute<3>()->info().is_empty );
          CGAL_assertion( !top_2->template attribute<3>()->info().is_empty );
          CGAL_assertion( !down_1->template attribute<3>()->info().is_empty );
          CGAL_assertion( !down_2->template attribute<3>()->info().is_empty );
          
          // poly_first  - poly_second            = {0}
          // poly_second - poly_first             = {0}
          // poly_first \cap poly_second          = P1P2 or Q1Q2 or P1Q1 or Q1P2
          // opposite( poly_first U poly_second ) = P2P1 or Q2Q1 or P2Q1 or Q2P1
          sew_2_marked_darts( final_map(),top_1    , top_2     ,mark_index, nodes, indices, polyline_info); //P1P2 or Q1Q2 or P1Q1 or Q1P2
          sew_2_marked_darts( final_map(),down_2   , down_1    ,mark_index, nodes, indices, polyline_info); //P2P1 or Q2Q1 or P2Q1 or Q2P1
          sew_3_marked_darts( final_map(),not_top_1, not_down_1 ,top_1,down_1,mark_index,darts_to_remove);
          sew_3_marked_darts( final_map(),not_top_2, not_down_2 ,top_2,down_2,mark_index,darts_to_remove);
          return true;
        }
      }
      sew_2_four_volumes_case_1(first_hedge->opposite(),second_hedge->opposite(),std::make_pair(indices.second,indices.first),nodes,index_p2,index_p1,index_q2,selected_hedge_to_dart,mark_index,darts_to_remove,polyline_info);
      return true;
    }
    if (index_p1==index_q2){
      if(index_p2!=-1){
        CGAL_assertion(index_p2!=index_q2);
        if(index_p2==index_q1){
          //-----first polyhedron
          Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_P2_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
          //-----second polyhedron
          Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);  
          
          Dart_handle top_1=dof_P1_outside->beta(3), not_top_1=dof_Q2_outside;
          Dart_handle top_2=dof_P2_outside->beta(3), not_top_2=dof_Q1_outside;
          Dart_handle down_1=dof_P1_outside, not_down_1=dof_Q2_outside->beta(3);
          Dart_handle down_2=dof_P2_outside, not_down_2=dof_Q1_outside->beta(3);
          if ( top_1->template attribute<3>()->info().is_empty ) std::swap(top_1,not_top_1);
          if ( top_2->template attribute<3>()->info().is_empty ) std::swap(top_2,not_top_2);
          if ( down_1->template attribute<3>()->info().is_empty ) std::swap(down_1,not_down_1);
          if ( down_2->template attribute<3>()->info().is_empty ) std::swap(down_2,not_down_2);
          CGAL_assertion( !top_1->template attribute<3>()->info().is_empty );
          CGAL_assertion( !top_2->template attribute<3>()->info().is_empty );
          CGAL_assertion( !down_1->template attribute<3>()->info().is_empty );
          CGAL_assertion( !down_2->template attribute<3>()->info().is_empty );          
          
          // poly_first  - poly_second            = P1P2 or P1Q1 or Q2P2 or Q2Q1
          // poly_second - poly_first             = Q1Q2 or Q1P1 or P2P1 or P2Q2
          // poly_first \cap poly_second          = {0}
          // opposite( poly_first U poly_second ) = all space
          sew_2_marked_darts( final_map(),top_1    , top_2     ,mark_index, nodes, indices, polyline_info); //P1P2 or Q1Q2 or P1Q1 or Q1P2
          sew_2_marked_darts( final_map(),down_2   , down_1    ,mark_index, nodes, indices, polyline_info); //P2P1 or Q2Q1 or P2Q1 or Q2P1
          sew_3_marked_darts( final_map(),not_top_1, not_down_1 ,top_1,down_1,mark_index,darts_to_remove);
          sew_3_marked_darts( final_map(),not_top_2, not_down_2 ,top_2,down_2,mark_index,darts_to_remove);
          return true;           
        }
      }
      sew_2_four_volumes_case_1(first_hedge->opposite(),second_hedge,std::make_pair(indices.second,indices.first),nodes,index_p2,index_p1,index_q1,selected_hedge_to_dart,mark_index,darts_to_remove,polyline_info,true);              
      return true;
    }
  }
  
  if(index_p2!=-1){
    if (index_p2==index_q1){
      sew_2_four_volumes_case_1(first_hedge,second_hedge->opposite(),indices,nodes,index_p1,index_p2,index_q2,selected_hedge_to_dart,mark_index,darts_to_remove,polyline_info,true);
      return true;
    }
    if(index_p2==index_q2){
      sew_2_four_volumes_case_1(first_hedge,second_hedge,indices,nodes,index_p1,index_p2,index_q1,selected_hedge_to_dart,mark_index,darts_to_remove,polyline_info);
      return true;
    }
  }
  
  
  return false;
}

//===//
//end//
//===//
  bool do_not_build_cmap; //set to true in the case only the corefinement must be done
  int number_coplanar_vertices; //number of intersection points between coplanar facets, see fixes XSL_TAG_CPL_VERT
public:
  Node_visitor_refine_polyhedra (Combinatorial_map_3_* ptr=NULL,bool do_not_build_cmap_=false):do_not_build_cmap(do_not_build_cmap_)
  {
    if (ptr!=NULL){
      final_map_comes_from_outside=true;
      final_map_ptr=ptr;
    }
    else{
      final_map_comes_from_outside=false;
      final_map_ptr=new Combinatorial_map_3_();
    }
  }
  
  ~Node_visitor_refine_polyhedra(){
    if(!final_map_comes_from_outside) delete final_map_ptr;
  }

  
  typedef internal_IOP::Predicates_on_constructions Node_storage_type;
  typedef Tag_false Is_polyhedron_const;
  static const bool do_need_vertex_graph = true;  //because we need to know which edges are constrained

  typedef Combinatorial_map_3_ Combinatorial_map_3;
  typedef internal_IOP::Volume_info<Polyhedron> Volume_info;

  const Combinatorial_map_3& combinatorial_map() const {return *final_map_ptr;}
  
  void set_number_of_intersection_points_from_coplanar_facets(int n){
    number_coplanar_vertices=n;
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
        typename In_halfedge_map::iterator it_hedge_map=in_hedge.insert(std::make_pair(additional_edge,Node_id_set())).first;
        it_hedge_map->second.insert(node_id);
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
        typename In_halfedge_map::iterator it_hedge_map=in_hedge.insert(std::make_pair(principal_edge,Node_id_set())).first;
        it_hedge_map->second.insert(node_id);
        check_node_on_non_manifold_edge(node_id,principal_edge);
      }
    }
  }
  
  template<class Iterator>
  void annotate_graph(Iterator begin,Iterator end)
  {
//    std::cout << "Annotation graph..." << std::endl;
    for (Iterator it=begin;it!=end;++it)
    {
      int node_id=it->first;
      if (non_manifold_nodes.find(node_id)!=non_manifold_nodes.end()) it->second.make_terminal();
      const std::set<int>& neighbors = it->second.neighbors;
      graph_of_constraints.insert(std::make_pair(node_id,neighbors));
    }
  }
  
  void update_terminal_nodes(std::vector<bool>&)
  {
    CGAL_assertion(!"Must not call this function");
  }
  
  void add_filtered_intersection(Halfedge_handle eh,Halfedge_handle fh,Polyhedron& Pe,Polyhedron& Pf){
    //use the representant halfedge of the facet as key
    //--set polyhedron for the two facets incident to the edge
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
    
    //used when object was create with hedge but opposite was used to split the original face
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
    #ifndef NDEBUG
    std::pair<typename Poly_to_map_node::iterator,bool> res = 
    #endif
      polyhedron_to_map_node_to_polyhedron_vertex.insert(std::make_pair( &P,Node_to_polyhedron_vertex_map() ));
    #ifndef NDEBUG
    CGAL_assertion(res.second == true);
    #endif    
  }
  
  //1) split_halfedges and retriangulate faces with no intersection point interior to the facet
  //2) retriangulate using a constrained Delaunay triangulation each triangle in each Polyhedron that contains at least
  //   one intersection point inside the facet
  //3) mark polyhedron edges that are on the intersection
  //4) create one output polyhedron per connected component of polyhedron, connected by an edge which is not an intersection edge   
  //5) import each piece into a common combinatorial map
  //6) glue all the pieces together
  template <class Node_vector>
  void finalize(const Node_vector& nodes){
    //mark halfedge that are on the intersection
    //SL: I needed to use a map because to get the orientation around the edge,
    //    I need to know in the case the third vertex is a node its index (for exact construction)
    typedef std::map<Halfedge_const_handle,std::pair<int,int>,Cmp_unik_ad > Border_halfedges_map;
    Border_halfedges_map border_halfedges;
    
    //Additionnal stuct to mark halfedge of the original polyhedra that are on the intersection
    typedef Halfedge_marker<Mark_intersection_halfedges> Marker;
    
    //store for each triangle facet which boundary is intersected by the other surface,
    //original vertices (and halfedges in the refined mesh pointing on these vertices)
    typedef std::map<Face_handle,Polyhedron_face_boundary,Cmp_handle> Faces_boundary;
    Faces_boundary faces_boundary; 
    
    //0) For each polyhedron, collect original vertices that belongs to the intersection.
    //   From the graph of constaints, extract intersection edges that are incident to such vertices. In case
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
              #ifndef NDEBUG
              Halfedge_handle start=hedge;
              #endif
              while ( hedge->opposite()->vertex()!=it_node_2_hedge_two->second->vertex() ){
                hedge=hedge->next()->opposite();
                #ifndef NDEBUG
                CGAL_assertion(hedge!=start);
                #endif
              }
              std::pair<int,int> edge_pair(*it_id,node_id_of_first);
              border_halfedges.insert( std::make_pair(hedge,edge_pair) );
              Marker::mark(hedge);
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
      Node_ids node_ids;   //indices of the intersection points to be inserted
      //we used a set to avoid having duplicated nodes reported on an edge of two coplanar triangles
      std::copy(it->second.begin(),it->second.end(),std::back_inserter(node_ids));
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
            
      #ifndef NDEBUG
      Vertex_handle original_vertex=hedge->opposite()->vertex();
      #endif
      
      //We need an edge incident to the source vertex of hedge. This is the first opposite edge created.      
      bool first=true; Halfedge_handle hedge_incident_to_src;
      //do split the edges
      for (std::vector<int>::const_iterator it_id=node_ids.begin();it_id!=node_ids.end();++it_id){
        Vertex_handle v=split_edge(hedge,nodes[*it_id],*P);
        node_to_polyhedron_vertex.insert(std::make_pair(*it_id,v));
        if (first){
          first=false;
          hedge_incident_to_src=hedge->opposite()->next();
        }
      }
      
      #ifndef NDEBUG
      CGAL_assertion(hedge_incident_to_src->vertex()==original_vertex);
      CGAL_assertion(hedge_incident_to_src->face()==hedge->opposite()->face());
      #endif

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
      typename Kernel::Plane_3 plane(triangle_boundary[0]->point(),triangle_boundary[1]->point(),triangle_boundary[2]->point());
      #else
      CGAL::Cartesian_converter<Kernel,Exact_kernel> convert;
      typename Exact_kernel::Plane_3 plane(convert(triangle_boundary[0]->point()),convert(triangle_boundary[1]->point()),convert(triangle_boundary[2]->point()));
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
      triangle_vertices[0]=triangulation.insert(plane.to_2d(triangle_boundary[0]->point()));
      triangle_vertices[1]=triangulation.insert(plane.to_2d(triangle_boundary[1]->point()));
      triangle_vertices[2]=triangulation.insert(plane.to_2d(triangle_boundary[2]->point()));
      #else
      //we can do this because these are input points.
      triangle_vertices[0]=triangulation.insert(plane.to_2d(convert(triangle_boundary[0]->point())));
      triangle_vertices[1]=triangulation.insert(plane.to_2d(convert(triangle_boundary[1]->point())));
      triangle_vertices[2]=triangulation.insert(plane.to_2d(convert(triangle_boundary[2]->point())));      
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
            res->second.insert( vn->info() ).second;
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
      internal_IOP::Triangulate_a_face<typename Polyhedron::HalfedgeDS> modifier(f,nodes,node_ids,node_to_polyhedron_vertex,edge_to_hedge,triangulation);
      
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
          border_halfedges.insert( std::make_pair(Halfedge_const_handle(it_poly_hedge->second),*it_cst) );
          Marker::mark(it_poly_hedge->second);
          update_edge_per_polyline(P,it_poly_hedge->first,it_poly_hedge->second);
        }
        else{
          //WARNING: in few case this is needed if the marked edge is on the border
          //to optimize it might be better to only use sorted pair. TAG_SLXX1
          std::pair<int,int> opposite_pair(it_cst->second,it_cst->first);
          it_poly_hedge=edge_to_hedge.find(opposite_pair);
          CGAL_assertion( it_poly_hedge!=edge_to_hedge.end() );

          border_halfedges.insert( std::make_pair(Halfedge_const_handle(it_poly_hedge->second),opposite_pair) );
          Marker::mark(it_poly_hedge->second);
          update_edge_per_polyline(P,it_poly_hedge->first,it_poly_hedge->second);          
        }
      }
    }
   
    if (do_not_build_cmap) return;
    
    //4) create one output polyhedron per connected component of polyhedron,
    //   connected by an edge which is not an intersection edge  
    //5) import into a Combinatorial map
    #ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "Nb marked edges " << border_halfedges.size() << std::endl;
//    for (typename Border_halfedges_map::iterator it=border_halfedges.begin();it!=border_halfedges.end();++it)
//      std::cout << it->first->opposite()->vertex()->point() << " " << it->first->vertex()->point() << " is constrained " << std::endl;
    std::cout << "Nb polylines " << an_edge_per_polyline.size() << std::endl;
    #endif
    
    internal_IOP::Non_intersection_halfedge<Polyhedron> criterium(border_halfedges);

    int mark_index=final_map().get_new_mark(); //mark used to tag dart that are on an intersection
    
    //define a map that will contain the correspondance between selected halfedges of the boundary and
    //their corresponding Dart_handle in the cmap.
    typedef std::map<Halfedge_const_handle,typename Combinatorial_map_3::Dart_handle,internal_IOP::Compare_address<Polyhedron> > Halfedge_to_dart_map;
    Halfedge_to_dart_map selected_hedge_to_dart;
    for (typename An_edge_per_polyline_map::iterator it=an_edge_per_polyline.begin();it!=an_edge_per_polyline.end();++it)
    {
      CGAL_assertion(it->second.first.size()==2);
      //orientation of faces around the edge (to be sure we can do it)
      Halfedge_handle first_hedge=it->second.first.begin()->second;
      Halfedge_handle second_hedge=boost::next(it->second.first.begin())->second;
      
      if (!first_hedge->is_border())               selected_hedge_to_dart.insert(std::make_pair(first_hedge,Dart_handle(NULL)));
      if (!first_hedge->opposite()->is_border())   selected_hedge_to_dart.insert(std::make_pair(first_hedge->opposite(),Dart_handle(NULL)));
      if (!second_hedge->is_border())              selected_hedge_to_dart.insert(std::make_pair(second_hedge,Dart_handle(NULL)));
      if (!second_hedge->opposite()->is_border())  selected_hedge_to_dart.insert(std::make_pair(second_hedge->opposite(),Dart_handle(NULL)));
    }
    
    #ifdef CGAL_COREFINEMENT_DEBUG
    int polynb=0;
    #endif
    for (typename Poly_to_map_node::iterator 
          it=polyhedron_to_map_node_to_polyhedron_vertex.begin();
          it!=polyhedron_to_map_node_to_polyhedron_vertex.end();
          ++it
        )
    {
      typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
      typedef ::CGAL::Union_find<Facet_const_handle> UF;
      typedef typename UF::handle UF_handle;
      typedef std::map<Facet_const_handle,std::list<Facet_const_handle>,internal::Compare_handle_ptr<Polyhedron> > Result;
      typedef std::map<Facet_const_handle,UF_handle,internal::Compare_handle_ptr<Polyhedron> > Facet_to_handle_map;
      
      UF uf;
      Facet_to_handle_map map_f2h;
      Result result;
      Polyhedron* current_poly=it->first;
      
      #ifdef CGAL_COREFINEMENT_DEBUG
      std::cout << "writing poly debug"<< std::endl;
      std::stringstream ss; 
      ss << "output_debug-" << ++polynb << ".off";
      std::ofstream output_debug(ss.str().c_str());
      output_debug << *current_poly;
      #endif
      
      extract_connected_components(*(current_poly),criterium,uf,map_f2h,result);

      
      //add each connected component in the map with 2 volumes per component.
      for (typename Result::iterator it_res=result.begin();it_res!=result.end();++it_res)
      {
        //create in the final Cmap a 2D component containing faces of a connected component 
        //(twice: one with same orientation and one with the opposite orientation to model the other volume)
        Dart_handle d = import_from_polyhedron_subset<Polyhedron>(final_map(),it_res->second.begin(),it_res->second.end(),criterium,selected_hedge_to_dart,mark_index);
        //set an attribute to one volume represented by this component to indicates 
        //a part outside of the polyhedron current_poly
        typename Combinatorial_map_3_::template Attribute_range<3>::type::iterator attrib=final_map().template create_attribute<3>();
        attrib->info().outside.insert(current_poly);
        final_map().template set_attribute<3>(d,attrib);
        //set the attribute for the opposite volume: represent a part inside current_poly
        attrib=final_map().template create_attribute<3>();
        attrib->info().inside.insert(current_poly);
        final_map().template set_attribute<3>(d->beta(3),attrib);        
        
        #ifdef CGAL_COREFINEMENT_DEBUG
        final_map().display_characteristics(std::cout);
        std::cout << std::endl;
        #endif
      }
    }
    #ifndef NDEBUG
    for(typename Halfedge_to_dart_map::iterator it=selected_hedge_to_dart.begin();it!=selected_hedge_to_dart.end();++it)
      CGAL_assertion(it->second!=Dart_handle(NULL));
    #endif
    
    CGAL_assertion(final_map().is_valid());

    std::set<Dart_handle> darts_to_remove;
    
    //6) Glue pieces together
    //   using one edge per intersection polyline, we merge the different volumes
    for (typename An_edge_per_polyline_map::iterator it=an_edge_per_polyline.begin();it!=an_edge_per_polyline.end();++it)
    {
      CGAL_assertion(it->second.first.size()==2);
      //orientation of faces around the edge (to be sure we can do it)
      std::pair<int,int> indices = it->first;
      const std::pair<bool,int>& polyline_info=it->second.second;

      //get the two halfedges incident to the edge [indices.first,indices.second]
      Halfedge_handle first_hedge=it->second.first.begin()->second;
      Halfedge_handle second_hedge=boost::next(it->second.first.begin())->second;

      CGAL_assertion(nodes[indices.second]==first_hedge->vertex()->point());
      CGAL_assertion(nodes[indices.first]==first_hedge->opposite()->vertex()->point());
      CGAL_assertion(nodes[indices.second]==second_hedge->vertex()->point());
      CGAL_assertion(nodes[indices.first]==second_hedge->opposite()->vertex()->point());
      
      Polyhedron* first_poly  = it->second.first.begin()->first;
      Polyhedron* second_poly = boost::next(it->second.first.begin())->first;
      
      //different handling depending on the number of incident triangles to the edge.
      //After sewing there are two,three or four volumes if there are two,three or four incident triangles respectively
      if ( first_hedge->is_border() || first_hedge->opposite()->is_border() ){
        if (second_hedge->is_border()  || second_hedge->opposite()->is_border())
          sew_2_two_volumes_case(first_hedge,second_hedge,border_halfedges,selected_hedge_to_dart,mark_index,darts_to_remove,nodes,indices,polyline_info);
        else
          sew_2_three_volumes_case(second_hedge, first_hedge,indices,nodes,border_halfedges,selected_hedge_to_dart,second_poly,first_poly,mark_index,darts_to_remove,polyline_info);
      }
      else
        if (second_hedge->is_border()  || second_hedge->opposite()->is_border())
          sew_2_three_volumes_case(first_hedge, second_hedge,indices,nodes,border_halfedges,selected_hedge_to_dart,first_poly,second_poly,mark_index,darts_to_remove,polyline_info);
        else
        {
          //Sort the four triangle facets around their common edge
          //  we suppose that the exterior of the polyhedron is indicated by
          //  counterclockwise oriented facets.
          Vertex_handle P1=first_hedge->opposite()->next()->vertex();
          Vertex_handle P2=first_hedge->next()->vertex();
          //    when looking from the side of indices.second, the interior of the first polyhedron is described 
          //    by turning counterclockwise from P1 to P2
          Vertex_handle Q1=second_hedge->opposite()->next()->vertex();
          Vertex_handle Q2=second_hedge->next()->vertex();
          //    when looking from the side of indices.second, the interior of the second polyhedron is described 
          //    by turning from Q1 to Q2
          
          //check if the third point of each triangular face is an original point (stay -1)
          //or a intersection point (in that case we need the index of the corresponding node to
          //have the exact value of the point)      
          int index_p1=node_index_of_incident_vertex(first_hedge->opposite()->next(),border_halfedges);
          int index_p2=node_index_of_incident_vertex(first_hedge->next(),border_halfedges);
          int index_q1=node_index_of_incident_vertex(second_hedge->opposite()->next(),border_halfedges);
          int index_q2=node_index_of_incident_vertex(second_hedge->next(),border_halfedges);      
          
          #ifdef CGAL_COREFINEMENT_DEBUG
          std::cout << index_p1 << " " << index_p2 << " " << index_q1 << " " <<index_q2 << std::endl;
          std::cout << nodes[indices.first] << " | " << nodes[indices.second] << std::endl;
          std::cout << P1->point() << " | " << P2->point() << " | " << Q1->point() << " | " <<Q2->point() << std::endl;
          #endif
          
          if ( coplanar_triangles_case_handled(first_hedge,second_hedge,indices,nodes,index_p1,index_p2,index_q1,index_q2,selected_hedge_to_dart,mark_index,darts_to_remove,polyline_info) )
            continue;
          
          CGAL_assertion(P1->point() !=Q1->point() && P1->point()!=Q2->point() && P2->point() !=Q1->point() && P2->point()!=Q2->point());
          
          bool Q1_is_between_P1P2 = filtered_order_around_edge(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes);
          bool Q2_is_between_P1P2 = filtered_order_around_edge(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes);
          
          
          //Recover the dart that will be the start point of the different sewing
          //  dof_X_outside = dart of face of , meaning the triangle containing the
          //  point X and part of the volume outside of the corresponding polyhedron      
          //-----first polyhedron
          Dart_handle dof_P1_outside = get_associated_dart(first_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_P2_outside = get_associated_dart(first_hedge,selected_hedge_to_dart);
          //-----second polyhedron
          Dart_handle dof_Q1_outside = get_associated_dart(second_hedge->opposite(),selected_hedge_to_dart);
          Dart_handle dof_Q2_outside = get_associated_dart(second_hedge,selected_hedge_to_dart);   
          
          if ( Q1_is_between_P1P2 ){
            if( Q2_is_between_P1P2 )
            {
              bool P1_is_between_Q1Q2 = filtered_order_around_edge(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes);
              if (!P1_is_between_Q1Q2){
                // poly_first  - poly_second            = P1Q1 U Q2P2
                // poly_second - poly_first             = {0}
                // poly_first \cap poly_second          = Q1Q2
                // opposite( poly_first U poly_second ) = P2P1
                sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //P1Q1
                sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q2P2
                sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q1Q2
                sew_2_marked_darts( final_map(),dof_P2_outside          , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //P2P1
                //update inside outside info (because darts from the same volume have been merged)
                dof_Q1_outside->beta(3)->template attribute<3>()->info().inside.insert(first_poly); //update Q1Q2 inside poly
                dof_P2_outside->template attribute<3>()->info().outside.insert(second_poly);//update P2P1 outside poly
              }
              else{
                // poly_first  - poly_second            = Q2Q1
                // poly_second - poly_first             = P2P1
                // poly_first \cap poly_second          = P1Q2 U Q1P2
                // opposite( poly_first U poly_second ) = {O}
                sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //Q2Q1
                sew_2_marked_darts( final_map(),dof_P2_outside          , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //P2P1
                sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q1P2
                sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P1Q2
                //update inside outside info (because darts from the same volume have been merged)
                dof_Q2_outside->template attribute<3>()->info().inside.insert(first_poly); //update Q2Q1 inside poly
                dof_P2_outside->template attribute<3>()->info().inside.insert(second_poly);//update P2P1 inside poly            
              }
            }
            else
            {
              // poly_first  - poly_second            = P1Q1
              // poly_second - poly_first             = P2Q2
              // poly_first \cap poly_second          = Q1P2
              // opposite( poly_first U poly_second ) = Q2P1
              sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //P1Q1
              sew_2_marked_darts( final_map(),dof_P2_outside          , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P2Q2
              sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q1P2
              sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q2P1
            }        
          }
          else
          {
            if( Q2_is_between_P1P2 )
            {
              // poly_first  - poly_second            = Q2P2
              // poly_second - poly_first             = Q1P1
              // poly_first \cap poly_second          = P1Q2
              // opposite( poly_first U poly_second ) = P2Q1
              sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q2P2
              sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q1P1
              sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P1Q2
              sew_2_marked_darts( final_map(),dof_P2_outside          , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //P2Q1
            }
            else
            {
              bool P1_is_between_Q1Q2 = filtered_order_around_edge(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes);
              if (!P1_is_between_Q1Q2){
                // poly_first  - poly_second            = P1P2
                // poly_second - poly_first             = Q1Q2
                // poly_first \cap poly_second          = {0}
                // opposite( poly_first U poly_second ) = P2Q1 U Q2P1
                sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P1P2
                sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //Q1Q2
                sew_2_marked_darts( final_map(),dof_P2_outside          , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //P2Q1
                sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q2P1
                //update inside outside info (because darts from the same volume have been merged)
                dof_Q1_outside->beta(3)->template attribute<3>()->info().outside.insert(first_poly); //update Q1Q2 outside poly
                dof_P1_outside->beta(3)->template attribute<3>()->info().outside.insert(second_poly);//update P2P1 outside poly            
              }
              else{
                // poly_first  - poly_second            = {0}
                // poly_second - poly_first             = Q1P1 U P2Q2
                // poly_first \cap poly_second          = P1P2
                // opposite( poly_first U poly_second ) = Q2Q1
                sew_2_marked_darts( final_map(),dof_Q1_outside->beta(3) , dof_P1_outside          ,mark_index, nodes, indices, polyline_info); //Q1P1
                sew_2_marked_darts( final_map(),dof_P2_outside          , dof_Q2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P2Q2
                sew_2_marked_darts( final_map(),dof_P1_outside->beta(3) , dof_P2_outside->beta(3) ,mark_index, nodes, indices, polyline_info); //P1P2
                sew_2_marked_darts( final_map(),dof_Q2_outside          , dof_Q1_outside          ,mark_index, nodes, indices, polyline_info); //Q2Q1
                //update inside outside info (because darts from the same volume have been merged)
                dof_P1_outside->beta(3)->template attribute<3>()->info().inside.insert(second_poly); //update P1P2 inside poly
                dof_Q2_outside->template attribute<3>()->info().outside.insert(first_poly);//update Q2Q1 outside poly            
              }
            }
          }
        }      
    }
    
    #ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "number of darts to remove: " << darts_to_remove.size() <<std::endl;
    #endif
    //remove darts from empty volumes
    for (typename std::set<Dart_handle>::iterator itdart=darts_to_remove.begin(),end=darts_to_remove.end();itdart!=end;++itdart){
      final_map().erase_dart(*itdart);
    }
    
    //remove empty volumes
    typedef typename Combinatorial_map_3_::template Attribute_range<3>::type Volume_attribute_range;
    Volume_attribute_range& ratrib=final_map().template attributes<3>();     
    typename Volume_attribute_range::iterator curr=ratrib.begin(),end=ratrib.end();
    do{
      if (curr->info().is_empty)
        final_map().template erase_attribute<3>(curr++);
      else
        ++curr;
    }
    while(curr!=end);
    
    CGAL_assertion(final_map().is_valid());
    
    //update the info of each volume knowing about only one polyhedron:
    //this happens when one polyhedron has a connected component
    //that do not intersect the other polyhedron
    typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel > Mesh_domain;
    CGAL_precondition(polyhedron_to_map_node_to_polyhedron_vertex.size()==2);
    Polyhedron* Poly_A = polyhedron_to_map_node_to_polyhedron_vertex.begin()->first;
    Polyhedron* Poly_B = boost::next(polyhedron_to_map_node_to_polyhedron_vertex.begin())->first;
    Mesh_domain* domain_A_ptr=NULL;
    Mesh_domain* domain_B_ptr=NULL;
    
    #ifdef CGAL_COREFINEMENT_DEBUG
    final_map().display_characteristics(std::cout); std::cout << "\n";
    #endif
    
    typename Combinatorial_map_3::template  One_dart_per_cell_range<3> cell_range=final_map().template one_dart_per_cell<3>();
    for (typename Combinatorial_map_3::template  One_dart_per_cell_range<3>::iterator
      it = cell_range.begin(), it_end=cell_range.end();
      it_end!= it;
      ++it )
    {
      internal_IOP::Volume_info<Polyhedron>& info=it->template attribute<3>()->info();
      std::size_t inside_size=info.inside.size();
      std::size_t outside_size=info.outside.size();
      if ( inside_size + outside_size == 1)
      {
        bool is_inside = (inside_size==1);
        Polyhedron* current_poly= is_inside? (*info.inside.begin()):(*info.outside.begin());
        Polyhedron* test_poly; 
        Mesh_domain* domain_ptr;
        if ( current_poly==Poly_A)
        {
          test_poly=Poly_B;
          if (domain_B_ptr == NULL) domain_B_ptr=new Mesh_domain(*Poly_B);
          domain_ptr=domain_B_ptr;
        }
        else
        {
          test_poly=Poly_A;
          if (domain_A_ptr == NULL) domain_A_ptr=new Mesh_domain(*Poly_A);
          domain_ptr=domain_A_ptr;
        }
        typename Mesh_domain::Is_in_domain is_in_domain(*domain_ptr);
        typename Kernel::Point_3 query=it->template attribute<0>()->point();
        if ( is_in_domain(query) )
          info.inside.insert(test_poly);
        else
          info.outside.insert(test_poly);
      }

      #ifdef CGAL_COREFINEMENT_DEBUG
      std::cout << "This volume has inside: ";
      for (typename std::set<Polyhedron*>::iterator itpoly=info.inside.begin();itpoly!=info.inside.end();++itpoly)
        std::cout << " " << *itpoly;
      std::cout << " and outside: ";
      for (typename std::set<Polyhedron*>::iterator itpoly=info.outside.begin();itpoly!=info.outside.end();++itpoly)
        std::cout << " " << *itpoly;
      std::cout << std::endl;
      #endif
    }
    if (domain_A_ptr!=NULL) delete domain_A_ptr;
    if (domain_B_ptr!=NULL) delete domain_B_ptr;
    
  }
};

}//namespace CGAL
  


#endif //CGAL_INTERSECTION_OF_POLYHEDRA_3_REFINEMENT_VISITOR_H
