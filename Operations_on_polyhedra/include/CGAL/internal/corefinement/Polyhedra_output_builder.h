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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_COREFINEMENT_POLYHEDRA_OUTPUT_BUILDER_H
#define CGAL_INTERNAL_COREFINEMENT_POLYHEDRA_OUTPUT_BUILDER_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/internal/corefinement/predicates.h>
#include <CGAL/internal/corefinement/connected_components.h>
#include <CGAL/internal/corefinement/utils.h>
#include <CGAL/array.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>

#include <bitset>

#include <sstream>
#include <fstream>

namespace CGAL{

namespace internal_IOP{

#ifdef DEFINE_UNUSED_CLASSES
template <class Polyhedron, bool reverse_patch_orientation>
struct Import_patch_helper;

template <class Polyhedron>
struct Import_patch_helper<Polyhedron, false>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;

  Vertex_handle first_vertex(Face_handle f){
    return f->halfedge()->vertex();
  }
  Vertex_handle second_hedge(Face_handle f){
    return f->halfedge()->next()->vertex();
  }
  Vertex_handle third_hedge(Face_handle f){
    return f->halfedge()->prev()->vertex();
  }
};

template <class Polyhedron>
struct Import_patch_helper<Polyhedron, true>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;

  Vertex_handle first_vertex(Face_handle f){
    return f->halfedge()->vertex();
  }
  Vertex_handle second_vertex(Face_handle f){
    return f->halfedge()->prev()->vertex();
  }
  Vertex_handle third_vertex(Face_handle f){
    return f->halfedge()->next()->vertex();
  }
};

template <class Polyhedron, bool reverse_patch_orientation = false>
class Import_patch
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  , public Import_patch_helper<Polyhedron, reverse_patch_orientation>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;

  const std::vector<Face_handle>& facets;
  const std::set<Vertex_handle>& interior_vertices;
  const std::vector<Halfedge_handle>& interior_halfedges;
  const std::vector<Halfedge_handle>& patch_border_halfedges;

public:

  Import_patch(
    const std::vector<Face_handle>& facets_,
    const std::set<Vertex_handle>& interior_vertices_,
    const std::vector<Halfedge_handle>& interior_halfedges_,
    const std::vector<Halfedge_handle>& patch_border_halfedges_
  ) : facets(facets_)
    , interior_vertices(interior_vertices_)
    ,interior_halfedges(interior_halfedges_)
    ,patch_border_halfedges(patch_border_halfedges_)
  {}

  void operator()(HDS& hds)
  {
    Polyhedron_incremental_builder_3<HDS> B( hds);
    B.begin_surface( interior_vertices.size(),
                     facets.size(),
                     interior_halfedges.size() * 2);

    std::size_t index=0;
    std::map<Vertex_handle,std::size_t> vertex_indices;

    //insert interior vertices
    BOOST_FOREACH( Vertex_handle vh, interior_vertices){
      B.add_vertex( vh->point() );
      vertex_indices.insert( std::make_pair(vh, index++) );
    }

    //insert border vertices : TODO WE SHOULD GET THEM FROM hds
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      if( vertex_indices.insert( std::make_pair(h->vertex(),index) ).second ){
        B.add_vertex( h->vertex()->point() );
        ++index;
      }
      if( vertex_indices.insert( std::make_pair(h->opposite()->vertex(),index) ).second ){
        B.add_vertex( h->opposite()->vertex()->point() );
        ++index;
      }
    }

    BOOST_FOREACH(Face_handle fh, facets)
    {
      B.begin_facet();
      B.add_vertex_to_facet( vertex_indices[ this->first_vertex(fh) ] );
      B.add_vertex_to_facet( vertex_indices[ this->second_vertex(fh) ] );
      B.add_vertex_to_facet( vertex_indices[ this->third_vertex(fh) ] );
      B.end_facet();
    }

    B.end_surface();
  }
};
#endif //DEFINE_UNUSED_CLASSES

template < class Polyhedron, class PolyhedronPointPMap, class HalfedgeOutputIterator >
struct Import_polyline
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;
  typedef std::map< Halfedge_handle,
                    Halfedge_handle,
                    Compare_unik_address<Polyhedron> > Hedge_map;
  typedef std::map< Vertex_handle,
                    Vertex_handle > Vertex_map;
  typedef internal_IOP::Compare_unik_address<Polyhedron>            Cmp_unik_ad;

  Vertex_map& vertex_map;
  Hedge_map&  P_to_O_hedge;
  Hedge_map&  Q_to_O_hedge;
  int nb_segments;
  Halfedge_handle P_first_halfedge;
  Halfedge_handle Q_first_halfedge;
  Halfedge_handle O_first_halfedge;
  const std::map< Halfedge_const_handle,
                  std::pair<int,int>,Cmp_unik_ad >& border_halfedges;
  PolyhedronPointPMap ppmap;
  HalfedgeOutputIterator output;

  Import_polyline(
    Vertex_map& vertex_map_,
    Hedge_map&  P_to_O_hedge_,
    Hedge_map&  Q_to_O_hedge_,
    int nb_segments_,
    Halfedge_handle P_first_halfedge_,
    Halfedge_handle Q_first_halfedge_,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges_,
    PolyhedronPointPMap ppmap_,
    HalfedgeOutputIterator output_)
  : vertex_map( vertex_map_ )
  , P_to_O_hedge( P_to_O_hedge_ )
  , Q_to_O_hedge( Q_to_O_hedge_ )
  , nb_segments( nb_segments_ )
  , P_first_halfedge( P_first_halfedge_ )
  , Q_first_halfedge( Q_first_halfedge_ )
  , border_halfedges( border_halfedges_ )
  , ppmap( ppmap_ )
  , output(output_)
  {}

  void operator()(HDS& hds)
  {
    HalfedgeDS_decorator<HDS> decorator(hds);

    typename HDS::Halfedge dummy_hedge; //dummy default constructed halfedge
    O_first_halfedge=hds.edges_push_back(dummy_hedge, dummy_hedge);
    *output++=O_first_halfedge;

    //make sure the first vertex does not already exist
    Vertex_handle source;
    std::pair< typename std::map< Vertex_handle, Vertex_handle >::iterator, bool > insert_res=
      vertex_map.insert( std::make_pair( P_first_halfedge->opposite()->vertex(), source ) );


    if( insert_res.second )
    {
      source = decorator.vertices_push_back( *(P_first_halfedge->opposite()->vertex()) );
      decorator.set_vertex_halfedge(source, O_first_halfedge->opposite());
      put(ppmap, source, get(ppmap, P_first_halfedge->opposite()->vertex() ) );
      insert_res.first->second = source;
    }
    else
      source = insert_res.first->second;

    //make sure the target vertex does not already exist if it is a polyline endpoint
    Vertex_handle target;
    if ( nb_segments==1 )
    {
      insert_res =
        vertex_map.insert( std::make_pair( P_first_halfedge->vertex(), target ) );
      if (insert_res.second)
      {
        target=decorator.vertices_push_back( *(P_first_halfedge->vertex()) );
        decorator.set_vertex_halfedge(target, O_first_halfedge);
        put(ppmap, target, get(ppmap, P_first_halfedge->vertex() ) );
        insert_res.first->second = target;
      }
      else
        target = insert_res.first->second;
    }
    else{
      target=decorator.vertices_push_back( *(P_first_halfedge->vertex()) );
      decorator.set_vertex_halfedge(target, O_first_halfedge);
      put(ppmap, target, get(ppmap, P_first_halfedge->vertex() ) );
    }


    //update source and target vertex of the edge created
    decorator.set_vertex(O_first_halfedge, target);
    decorator.set_vertex(O_first_halfedge->opposite(), source);


    Halfedge_handle O_previous=O_first_halfedge;
    Halfedge_handle P_previous=P_first_halfedge;
    Halfedge_handle Q_previous=Q_first_halfedge;

    //set the correspondance
    P_to_O_hedge.insert( std::make_pair(P_previous, O_previous) );
    Q_to_O_hedge.insert( std::make_pair(Q_previous, O_previous) );

    source=target;
    for (int i=1; i<nb_segments; ++i)
    {
      //create new edge
      Halfedge_handle O_hedge = hds.edges_push_back(dummy_hedge, dummy_hedge);
      *output++=O_hedge;
      //get the new edge
      Halfedge_handle P_hedge = next_marked_halfedge_around_target_vertex(P_previous, border_halfedges);
      Halfedge_handle Q_hedge = next_marked_halfedge_around_target_vertex(Q_previous, border_halfedges);
      //if this is the final segment, only create a target vertex if it does not exist
      if (i!=nb_segments-1)
      {
        target=decorator.vertices_push_back( *(P_hedge->vertex()) );
        decorator.set_vertex_halfedge(target, O_hedge);
        put(ppmap, target, get(ppmap, P_hedge->vertex() ) );
      }
      else{
        std::pair< typename std::map< Vertex_handle, Vertex_handle >::iterator, bool > insert_res=
          vertex_map.insert( std::make_pair(P_hedge->vertex(), Vertex_handle()) );
        if (insert_res.second)
        {
          target=decorator.vertices_push_back( *(P_hedge->vertex()) );
          decorator.set_vertex_halfedge(target, O_hedge);
          put(ppmap, target, get(ppmap, P_hedge->vertex() ) );
          insert_res.first->second = target;
        }
        else
          target = insert_res.first->second;
      }

      decorator.set_vertex(O_hedge, target);
      decorator.set_vertex(O_hedge->opposite(), source);

      O_previous=O_hedge;
      P_previous = P_hedge;
      Q_previous = Q_hedge;
      source = target;

      P_to_O_hedge.insert( std::make_pair(P_previous, O_previous) );
      Q_to_O_hedge.insert( std::make_pair(Q_previous, O_previous) );
    }
  }
};

template <class Polyhedron, bool reverse_patch_orientation>
struct Surface_extension_helper;

template <class Polyhedron>
struct Surface_extension_helper<Polyhedron, true>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef std::map< Halfedge_handle,
                    Halfedge_handle,
                    Compare_unik_address<Polyhedron> > Hedge_map;
  Hedge_map& Qhedge_to_Phedge;

  Surface_extension_helper(Hedge_map& qh_to_ph):Qhedge_to_Phedge(qh_to_ph){}

  Halfedge_handle get_hedge(Halfedge_handle qhedge)
  {
    CGAL_assertion( Qhedge_to_Phedge.find(qhedge)!=Qhedge_to_Phedge.end() );
    std::pair<Halfedge_handle, Halfedge_handle> key_and_value =
      *Qhedge_to_Phedge.find(qhedge);
    return key_and_value.first != qhedge
           ? key_and_value.second
           : key_and_value.second->opposite();
  }

  Halfedge_handle first_hedge(Face_handle f){
    return get_hedge( f->halfedge() );
  }
  Halfedge_handle second_hedge(Face_handle f){
    return get_hedge( f->halfedge()->prev() );
  }
  Halfedge_handle third_hedge(Face_handle f){
    return get_hedge( f->halfedge()->next() );
  }
};

template <class Polyhedron>
struct Surface_extension_helper<Polyhedron, false>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef std::map< Halfedge_handle,
                    Halfedge_handle,
                    Compare_unik_address<Polyhedron> > Hedge_map;
  Hedge_map& Qhedge_to_Phedge;

  Surface_extension_helper(Hedge_map& qh_to_ph):Qhedge_to_Phedge(qh_to_ph){}

  Halfedge_handle get_hedge(Halfedge_handle qhedge)
  {
    CGAL_assertion( Qhedge_to_Phedge.count(qhedge)==1 );
    std::pair<Halfedge_handle, Halfedge_handle> key_and_value =
      *Qhedge_to_Phedge.find(qhedge);
    return key_and_value.first == qhedge
           ? key_and_value.second
           : key_and_value.second->opposite();
  }

  Halfedge_handle first_hedge(Face_handle f){
    return get_hedge( f->halfedge() );
  }
  Halfedge_handle second_hedge(Face_handle f){
    return get_hedge( f->halfedge()->next() );
  }
  Halfedge_handle third_hedge(Face_handle f){
    return get_hedge( f->halfedge()->prev() );
  }
};



template <class Polyhedron, class PolyhedronPointPMap, bool reverse_patch_orientation=false>
class Surface_extension_by_patch_appending
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
  , public Surface_extension_helper<Polyhedron, reverse_patch_orientation>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;

  const std::vector<Face_handle>& facets;
  const std::vector<Halfedge_handle>& interior_halfedges;
  const std::vector<Halfedge_handle>& patch_border_halfedges;
  const std::set<Vertex_handle>& interior_vertices;
  PolyhedronPointPMap ppmap;

  typedef Surface_extension_helper<Polyhedron, reverse_patch_orientation> Base;
  using Base::first_hedge;
  using Base::second_hedge;
  using Base::third_hedge;
  using Base::get_hedge;
  using Base::Qhedge_to_Phedge;

public:

  Surface_extension_by_patch_appending(
    const std::vector<Face_handle>& facets_,
    const std::vector<Halfedge_handle>& interior_halfedges_,
    const std::vector<Halfedge_handle>& patch_border_halfedges_,
    const std::set<Vertex_handle>& interior_vertices_,
    typename Base::Hedge_map& Qhedge_to_Phedge_,
    PolyhedronPointPMap ppmap_
  ) : Base( Qhedge_to_Phedge_ )
    ,facets(facets_)
    ,interior_halfedges(interior_halfedges_)
    ,patch_border_halfedges(patch_border_halfedges_)
    ,interior_vertices(interior_vertices_)
    ,ppmap(ppmap_)
  {}

  void operator()(HDS& hds)
  {
    HalfedgeDS_decorator<HDS> decorator(hds);
    std::vector<Halfedge_handle> interior_vertex_halfedges;

    //insert interior halfedges and create interior vertices
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges){
      Halfedge_handle new_h = hds.edges_push_back( *h );
      Qhedge_to_Phedge[ h ] = new_h;
      // put new halfedges on the border of the mesh
      decorator.set_face(new_h, Face_handle());
      decorator.set_face(new_h->opposite(), Face_handle());

      //create a copy of interior vertices only once
      if (  h->vertex()->halfedge()==h && interior_vertices.count(h->vertex()) )
      {
        Vertex_handle v = decorator.vertices_push_back( *(h->vertex()) );
        decorator.set_vertex_halfedge(v, new_h);
        decorator.set_vertex(new_h, v);
        put(ppmap, v, get(ppmap, h->vertex() ) );
        interior_vertex_halfedges.push_back( new_h );
      }
      if (  h->opposite()->vertex()->halfedge()==h->opposite() &&
            interior_vertices.count(h->opposite()->vertex()) )
      {
        Vertex_handle v = decorator.vertices_push_back( *(h->opposite()->vertex()) );
        decorator.set_vertex_halfedge(v, new_h->opposite());
        decorator.set_vertex(new_h->opposite(), v);
        put(ppmap, v, get(ppmap, h->opposite()->vertex() ) );
        interior_vertex_halfedges.push_back( new_h->opposite() );
      }
    }

    //create facets and connect halfedges
    BOOST_FOREACH(Face_handle f, facets)
    {
      Halfedge_handle hedges[]={
        first_hedge(f), second_hedge(f), third_hedge(f)
      };

      Face_handle new_f = decorator.faces_push_back( *f );
      decorator.set_face_halfedge(new_f, hedges[0]);

      for (int i=0;i<3;++i)
      {
        hedges[i]->HBase::set_next(hedges[(i+1)%3]);
        decorator.set_prev(hedges[(i+1)%3], hedges[i]);
        decorator.set_face(hedges[i], new_f);
      }
    }

    // handle interior edges that are on the border of the mesh:
    // they do not have a prev/next pointer set since only the pointers
    // of patch interior halfedges part a face have been. In the following
    // (i) we set the next/prev pointer around interior vertices on the mesh
    // boundary and (ii) we collect interior mesh border halfedges incident to
    // a patch border vertex and set their next/prev pointer (possibly of
    // another patch)

    // Containers used for step (ii) for collecting mesh border halfedges
    // with source/target on an intersection polyline that needs it prev/next
    // pointer to be set
    std::vector<Halfedge_handle> border_halfedges_source_to_link;
    std::vector<Halfedge_handle> border_halfedges_target_to_link;
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
      if (h->is_border_edge())
      {
        if (!h->is_border()) h=h->opposite();

        Vertex_handle src = h->opposite()->vertex();
        Vertex_handle tgt = h->vertex();
        if (reverse_patch_orientation) std::swap(src, tgt);

        if ( !interior_vertices.count(src) )
          border_halfedges_source_to_link.push_back(get_hedge(h));
        if ( !interior_vertices.count(tgt) ){
          border_halfedges_target_to_link.push_back(get_hedge(h));
          continue; // since the next halfedge should not be in the same patch
        }
        CGAL_assertion( h->is_border() &&
                        h->prev()->is_border() && h->next()->is_border() );
        // step (i)
        Halfedge_handle qhedge=get_hedge(h);
        Halfedge_handle qnext = reverse_patch_orientation ?
                               get_hedge(h->prev()) : get_hedge(h->next());
        CGAL_assertion( qhedge->is_border() && qnext->is_border()  );
        qhedge->HBase::set_next(qnext);
        decorator.set_prev(qnext, qhedge);
      }
    // now the step (ii) we look for the candidate halfedge by turning around
    // the vertex in the direction of the interior of the patch
    BOOST_FOREACH(Halfedge_handle h, border_halfedges_target_to_link)
    {
      Halfedge_handle candidate=h->opposite()->prev()->opposite();
      CGAL_assertion_code(Halfedge_handle start=candidate);
      while (!candidate->is_border()){
        candidate=candidate->prev()->opposite();
        CGAL_assertion(candidate!=start);
      }
      h->HBase::set_next(candidate);
      decorator.set_prev(candidate,h);
    }
    BOOST_FOREACH(Halfedge_handle h, border_halfedges_source_to_link)
    {
      Halfedge_handle candidate=h->opposite()->next()->opposite();
      while (!candidate->is_border())
        candidate=candidate->next()->opposite();
      candidate->HBase::set_next(h);
      decorator.set_prev(h,candidate);
    }

    // For all interior vertices, update the vertex pointer
    // of all but the vertex halfedge
    BOOST_FOREACH(Halfedge_handle h, interior_vertex_halfedges)
    {
      Vertex_handle v=h->vertex();
      Halfedge_handle next_around_vertex=h;
      do{
        CGAL_assertion (next_around_vertex->next() != Halfedge_handle());
        next_around_vertex=next_around_vertex->next()->opposite();
        decorator.set_vertex(next_around_vertex, v);
      }while(h != next_around_vertex);
    }

    // For all patch boundary vertices, update the vertex pointer
    // of all but the vertex halfedge
    BOOST_FOREACH(Halfedge_handle qhedge, patch_border_halfedges)
    {
      //check for a halfedge pointing inside an already imported patch
      Halfedge_handle h = get_hedge(qhedge);
      CGAL_assertion( h->next()!=Halfedge_handle() );
      // update the pointers on the target
      Halfedge_handle next_around_target=h;
      Vertex_handle v=h->vertex();
      do{
        next_around_target = next_around_target->next()->opposite();
        decorator.set_vertex(next_around_target, v);
      }while( next_around_target->next()!=Halfedge_handle() &&
              next_around_target!=h && !next_around_target->is_border());
      // update the pointers on the source
      Halfedge_handle next_around_source=h->prev();
      CGAL_assertion(next_around_source!=Halfedge_handle());
      v = h->opposite()->vertex();
      do{
        decorator.set_vertex(next_around_source, v);
        next_around_source = next_around_source->opposite()->prev();
      }while( next_around_source!=Halfedge_handle() &&
              next_around_source!=h->opposite() &&
              !next_around_source->is_border());
    }
  }
};

template <class Polyhedron>
struct Intersection_polylines{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  const std::vector<Halfedge_handle>& P;
  const std::vector<Halfedge_handle>& Q;
  const std::vector<int>& lengths;
  boost::dynamic_bitset<> to_skip;
  boost::dynamic_bitset<> to_skip_in_P;
  boost::dynamic_bitset<> to_skip_in_Q;
  Intersection_polylines(
    const std::vector<Halfedge_handle>& P_polylines,
    const std::vector<Halfedge_handle>& Q_polylines,
    const std::vector<int>& lengths_
  ) : P( P_polylines )
    , Q( Q_polylines )
    , lengths( lengths_ )
    , to_skip(P.size(),false)
    , to_skip_in_P(P.size(),false)
    , to_skip_in_Q(P.size(),false)
  {}
};

template <class Polyhedron>
class Remove_patch_simplices : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;

  const std::vector<Face_handle>& facets;
  const std::set<Vertex_handle>& interior_vertices;
  const std::vector<Halfedge_handle>& interior_halfedges;
  const std::vector<Halfedge_handle>& patch_border_halfedges;

public:

  Remove_patch_simplices(
    const std::vector<Face_handle>& facets_,
    const std::set<Vertex_handle>& interior_vertices_,
    const std::vector<Halfedge_handle>& interior_halfedges_,
    const std::vector<Halfedge_handle>& patch_border_halfedges_
  ) : facets(facets_)
    , interior_vertices(interior_vertices_)
    ,interior_halfedges(interior_halfedges_)
    ,patch_border_halfedges(patch_border_halfedges_)
  {}

  void operator()(HDS& hds)
  {
    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);

    // put the halfedges on the boundary of the patch on the boundary of the polyhedron
    Face_handle border;
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
      decorator.set_face(h, border);

    // set next/prev relationship of border halfedges
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      Halfedge_handle next=h->next();
      while(!next->is_border())
        next=next->opposite()->next();
      h->HBase::set_next(next);
      decorator.set_prev(next,h);
      decorator.set_vertex_halfedge(h->vertex(),h);
    }

    // In case a ccb of the patch is not a cycle (the source and target vertices
    // are border vertices), the first halfedge of that ccb will not have its
    // prev pointer set correctly. To fix that, we consider all interior edges
    // and check for one that is on the border of the patch and that is incident
    // to a border vertex and use it to get the missing prev pointer.
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
      if(h->is_border_edge())
      {
        if (h->is_border()) h=h->opposite();
        if ( !interior_vertices.count(h->vertex()) )
        {
          // look for the halfedge belonging to patch_border_halfedge
          // having the prev pointer not correctly set
          Halfedge_handle next=h->next();
          while(!next->is_border())
            next=next->opposite()->next();
          CGAL_assertion( next->is_border() );//we marked it above!
          // now update the prev pointer
          Halfedge_handle prev=h->opposite()->prev();
          prev->HBase::set_next(h->next());
          decorator.set_prev(h->next(),prev);
          decorator.set_vertex_halfedge(prev->vertex(),prev);
        }
      }

    //now remove the simplices
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
      hds.edges_erase(h);
    BOOST_FOREACH(Face_handle f, facets)
      hds.faces_erase(f);
    BOOST_FOREACH(Vertex_handle v, interior_vertices)
      hds.vertices_erase(v);
  }
};

template <class Polyhedron>
class Remove_isolated_patch_simplices
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;

  const std::vector<Face_handle>& facets;
  const std::set<Vertex_handle>& interior_vertices;
  const std::vector<Halfedge_handle>& interior_halfedges;
  const std::vector<Halfedge_handle>& border_halfedges;

public:

  Remove_isolated_patch_simplices(
    const std::vector<Face_handle>& facets_,
    const std::set<Vertex_handle>& interior_vertices_,
    const std::vector<Halfedge_handle>& interior_halfedges_,
    const std::vector<Halfedge_handle>& border_halfedges_
  ) : facets(facets_)
    , interior_vertices(interior_vertices_)
    , interior_halfedges(interior_halfedges_)
    , border_halfedges(border_halfedges_)
  {}

  void operator()(HDS& hds)
  {
    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);
    //remove the simplices
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
      hds.edges_erase(h);
    // There is no shared halfedge between duplicated patches even
    // if they were before the duplication. Thus the erase that follows is safe.
    // However remember that vertices were not duplicated which is why their
    // removal is not handled here (still in use or to be removed in
    // remove_unused_polylines())
    BOOST_FOREACH(Halfedge_handle h, border_halfedges)
      hds.edges_erase(h);
    BOOST_FOREACH(Face_handle f, facets)
      hds.faces_erase(f);
    BOOST_FOREACH(Vertex_handle v, interior_vertices)
      hds.vertices_erase(v);
  }
};

template <class Polyhedron>
class Disconnect_patches : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename HDS::Halfedge_handle   Halfedge_handle;
  typedef typename HDS::Vertex_handle     Vertex_handle;
  typedef typename HDS::Face_handle       Face_handle;
  typedef typename HDS::Vertex            Vertex;
  typedef typename HDS::Halfedge          Halfedge;
  typedef typename HDS::Face              Face;
  typedef typename HDS::Halfedge::Base    HBase;

  const std::vector<Face_handle>& facets;
  const std::set<Vertex_handle>& interior_vertices;
  const std::vector<Halfedge_handle>& interior_halfedges;
  const std::vector<Halfedge_handle>& patch_border_halfedges;
        std::vector<Halfedge_handle>& new_patch_border;

public:

  Disconnect_patches(
    const std::vector<Face_handle>& facets_,
    const std::set<Vertex_handle>& interior_vertices_,
    const std::vector<Halfedge_handle>& interior_halfedges_,
    const std::vector<Halfedge_handle>& patch_border_halfedges_,
    std::vector<Halfedge_handle>& new_patch_border_
  ) : facets(facets_)
    , interior_vertices(interior_vertices_)
    , interior_halfedges(interior_halfedges_)
    , patch_border_halfedges(patch_border_halfedges_)
    , new_patch_border(new_patch_border_)
  {}

  void operator()(HDS& hds)
  {
    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);

    new_patch_border.reserve( patch_border_halfedges.size() );

    std::map<Halfedge_handle, Halfedge_handle> old_to_new;

    // put the halfedges on the boundary of the patch on the boundary of the polyhedron
    Face_handle border;
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      Halfedge_handle new_hedge = hds.edges_push_back(*h);
      new_patch_border.push_back(new_hedge);
      decorator.set_face(h, border);
      decorator.set_face(new_hedge->opposite(), border);
      old_to_new.insert( std::make_pair(h, new_hedge) );
    }

    // update next/prev pointer of new hedges in case it is one of the new hedges
    BOOST_FOREACH(Halfedge_handle h, new_patch_border)
    {
      if ( h->next()->is_border() ){
        h->HBase::set_next( old_to_new[h->next()] );
        decorator.set_prev( h->next(), h);
      }
    }

    // set next/prev pointers on the border of the neighbor patch
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      Halfedge_handle next=h->next();
      // check if not already done
      if ( !next->is_border() ){
        do{
          next=next->opposite()->next();
        } while(!next->is_border());
        h->HBase::set_next(next);
        decorator.set_prev(next,h);
      }

      // setting prev is only needed in case the polyhedron has a boundary
      // and the intersection polyline intersects its boundary
      if ( !h->prev()->is_border() ){
        Halfedge_handle prev=h->prev();
        do{
          prev=prev->opposite()->prev();
        } while( !prev->is_border() );
        prev->HBase::set_next(h);
        decorator.set_prev(h, prev);
      }

      CGAL_assertion( h->prev()->is_border() );

      decorator.set_vertex_halfedge(h->vertex(),h);
      decorator.set_vertex_halfedge(h->opposite()->vertex(),h->opposite()); //  only needed if the polyhedra is open
    }

    //update next/prev relationship inside the patch
    //to have a correct connectivity, and update facet halfedge pointer
    BOOST_FOREACH(Halfedge_handle h, new_patch_border)
    {
      if ( h->prev()->next() != h )
        h->prev()->HBase::set_next( h );
      if ( h->next()->prev() != h )
        decorator.set_prev(h->next(), h);
      decorator.set_face_halfedge(h->facet(), h);
    }

    // update next/prev pointers on the border of the patch
    BOOST_FOREACH(Halfedge_handle h, new_patch_border)
    {
      h=h->opposite();
      //set next pointer if not already set
      if ( h->next()->prev()!=h )
      {
        // we visit facets inside the patch we consider
        Halfedge_handle candidate = h->opposite()->prev()->opposite();
        while ( !candidate->is_border() )
          candidate = candidate->prev()->opposite();
        h->HBase::set_next(candidate);
        decorator.set_prev(candidate,h);
      }
      CGAL_assertion( h->next()->prev()== h );

      // set prev pointer if not already set
      if ( h->prev()->next() != h )
      {
        Halfedge_handle candidate = h->opposite()->next()->opposite();
        while ( !candidate->is_border() )
          candidate = candidate->next()->opposite();
        decorator.set_prev(h,candidate);
        candidate->HBase::set_next(h);
      }

      CGAL_assertion( h->next()->prev()== h );
      CGAL_assertion( h->prev()->is_border() );
      CGAL_assertion( h->next()->is_border() );
    }
  }
};


template<class Polyhedron>
struct Patch_description{
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  std::vector<Facet_handle> facets;
  std::set<Vertex_handle> interior_vertices;
  std::vector<Halfedge_handle> interior_halfedges;
  std::vector<Halfedge_handle> patch_border_halfedges;
  bool is_initialized;

  Patch_description(): is_initialized(false) {};
};

template <class Polyhedron, class Facet_index_pmap, class Is_marked_edge_map>
void extract_patch_simplices(
  std::size_t patch_id,
  Polyhedron& P,
  Facet_index_pmap facet_indices,
  const std::vector<std::size_t>& patch_ids,
  std::vector<typename Polyhedron::Facet_handle>& facets,
  std::set<typename Polyhedron::Vertex_handle>& interior_vertices,
  std::vector<typename Polyhedron::Halfedge_handle>& interior_halfedges,
  std::vector<typename Polyhedron::Halfedge_handle>& patch_border_halfedges,
  const Is_marked_edge_map& is_marked_edge
)
{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;

  for (typename Polyhedron::Facet_iterator fit=P.facets_begin(),
                                           fit_end=P.facets_end();
                                           fit!=fit_end; ++fit)
  {
    if ( patch_ids[ get(facet_indices, fit) ]==patch_id )
    {
      facets.push_back( fit );
      Halfedge_handle hedges[]={fit->halfedge(), fit->halfedge()->next(), fit->halfedge()->prev()};
      for (int i=0;i<3;++i)
      {
        if ( !is_marked_edge.count(hedges[i]) )
        {
          if ( hedges[i] < hedges[i]->opposite()
               || hedges[i]->opposite()->is_border())
            interior_halfedges.push_back( hedges[i] );
        }
        else
          patch_border_halfedges.push_back(hedges[i]);
      }
    }
  }

  std::set<Vertex_handle> border_vertices;
  for (std::size_t k=0, end=patch_border_halfedges.size(); k!=end; ++k)
  {
    border_vertices.insert( patch_border_halfedges[k]->vertex() );
    // if the model is not closed i.e. patch_border_halfedge is not cycle only
    border_vertices.insert( patch_border_halfedges[k]->opposite()->vertex() );
  }

  for (std::size_t k=0, end=interior_halfedges.size(); k!=end; ++k)
  {
    if ( border_vertices.find( interior_halfedges[k]->vertex() )
         == border_vertices.end() )
    {
      interior_vertices.insert( interior_halfedges[k]->vertex() );
    }
    if ( border_vertices.find( interior_halfedges[k]->opposite()->vertex() )
         == border_vertices.end() )
    {
      interior_vertices.insert( interior_halfedges[k]->opposite()->vertex() );
    }
  }
}

template <class Polyhedron, class Facet_index_pmap>
struct Patch_container{
//typedefs
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef internal_IOP::Compare_unik_address<Polyhedron>            Cmp_unik_ad;
// data members
  std::vector< Patch_description<Polyhedron> > patches;
// external data members
  Polyhedron* poly_ptr;
  const std::vector<std::size_t>& patch_ids;
  Facet_index_pmap facet_id_pmap;
  const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges;
// constructor
  Patch_container(
    Polyhedron* poly_ptr_,
    const std::vector<std::size_t>& patch_ids_,
    Facet_index_pmap facet_id_pmap_,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges_,
    std::size_t nb_patches
  ) : patches(nb_patches)
    , poly_ptr(poly_ptr_)
    , patch_ids(patch_ids_)
    , facet_id_pmap(facet_id_pmap_)
    , border_halfedges(border_halfedges_)
  {}

  Patch_description<Polyhedron>& operator[](std::size_t i) {
    if ( !patches[i].is_initialized )
    {
      extract_patch_simplices(
        i, *poly_ptr,
        facet_id_pmap, patch_ids,
        patches[i].facets, patches[i].interior_vertices,
        patches[i].interior_halfedges, patches[i].patch_border_halfedges,
        border_halfedges
      );

      patches[i].is_initialized=true;
    }
    return patches[i];
  }

  /// debug
  std::ostream& dump_patch(std::size_t i, std::ostream& out)
  {
    Patch_description<Polyhedron>& patch=this->operator[](i);
    out << "OFF\n" << patch.interior_vertices.size()+patch.patch_border_halfedges.size();
    out << " " << patch.facets.size() << " 0\n";
    std::map<typename Polyhedron::Vertex_handle, int> vertexid;
    int id=0;
    BOOST_FOREACH(typename Polyhedron::Vertex_handle vh, patch.interior_vertices)
    {
      vertexid[vh]=id++;
      out << vh->point() << "\n";
    }

    BOOST_FOREACH(typename Polyhedron::Halfedge_handle hh, patch.patch_border_halfedges)
    {
      vertexid[hh->vertex()]=id++;
      out << hh->vertex()->point() << "\n";
    }

    BOOST_FOREACH(typename Polyhedron::Facet_handle fh, patch.facets)
    {
      out << "3 " << vertexid[fh->halfedge()->vertex()] <<
             " "  << vertexid[fh->halfedge()->next()->vertex()] <<
             " "  << vertexid[fh->halfedge()->next()->next()->vertex()] << "\n";
    }

    return out;
  }

  void dump_patches(const boost::dynamic_bitset<>& selection, std::string prefix)
  {
    for (std::size_t i=selection.find_first();
                     i < selection.npos; i = selection.find_next(i))
    {
      std::stringstream ss;
      ss << prefix << "-" << i << ".off";
      std::ofstream output(ss.str().c_str());
      dump_patch(i, output);
    }
  }

};

} //end of namespace internal_IOP

template<class Polyhedron>
struct Default_facet_id_pmap
{
  typedef boost::read_write_property_map_tag                             category;
  typedef std::size_t                                                  value_type;
  typedef std::size_t&                                                 reference;
  typedef typename Polyhedron::Facet_const_handle key_type;

  friend std::size_t get(Default_facet_id_pmap, key_type f) { return f->id(); }
  friend void put(Default_facet_id_pmap, key_type f, std::size_t i) { const_cast<typename Polyhedron::Facet&>(*f).id()=i; }
};

namespace Corefinement
{

template <class Polyhedron,
          class Facet_id_pmap_=Default,
          class Kernel_=Default,
          class PolyhedronPointPMap_=Default,
          class EdgeMarkPropertyMap_=Default >
class Polyhedra_output_builder
{
//Default typedefs
  typedef typename Default::Get<
    PolyhedronPointPMap_,
    Default_polyhedron_ppmap<Polyhedron> >::type            PolyhedronPointPMap;
  typedef typename Default::Get<
    Kernel_,
    typename Kernel_traits<
      typename boost::property_traits<PolyhedronPointPMap>::value_type
    >::Kernel >::type                                                    Kernel;
  typedef typename Default::Get<
    Facet_id_pmap_,
    Default_facet_id_pmap<Polyhedron> >::type                     Facet_id_pmap;
  typedef typename Default::Get<EdgeMarkPropertyMap_,
    Corefinement::Dummy_edge_mark_property_map<Polyhedron> >::type
                                                            EdgeMarkPropertyMap;


public:
//Boolean operation indices
  enum Boolean_operation {P_UNION_Q = 0, P_INTER_Q, P_MINUS_Q, Q_MINUS_P, NONE };
private:
//Data members
  Polyhedron *P_ptr, *Q_ptr;
  cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output;
  Facet_id_pmap P_facet_id_pmap, Q_facet_id_pmap;
  PolyhedronPointPMap ppmap;
  EdgeMarkPropertyMap edge_mark_pmap;
  bool input_with_coplanar_facets;
  // bitset containing information about operations that cannot be
  // performed because of non-manifoldness or that is ambiguous
  // 0 = P+Q
  // 1 = P inter Q
  // 2 = P - Q
  // 3 = Q  - P
  std::bitset<4> impossible_operation;

  //Orientation of polyhedra
  bool is_P_inside_out;
  bool is_Q_inside_out;

//Polyhedron typedefs
  typedef typename Polyhedron::Halfedge_const_handle      Halfedge_const_handle;
  typedef typename Polyhedron::Halfedge_handle                  Halfedge_handle;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
  typedef typename Polyhedron::Vertex_const_handle          Vertex_const_handle;
  typedef typename Polyhedron::Facet_handle                        Facet_handle;
//Other typedefs
  typedef internal_IOP::Compare_unik_address<Polyhedron>            Cmp_unik_ad;
  typedef std::map< std::pair<int,int>,
                    std::pair< std::map<Polyhedron*,Halfedge_handle>,
                               std::pair<bool,int> > > An_edge_per_polyline_map;
  typedef std::map<int,Vertex_handle> Node_to_polyhedron_vertex_map;
  typedef std::map<Polyhedron*, Node_to_polyhedron_vertex_map> Poly_to_map_node;
  typedef internal_IOP::Intersection_polylines<Polyhedron> Intersection_polylines;
  typedef internal_IOP::Patch_description<Polyhedron> Patch_description;
  typedef internal_IOP::Patch_container<Polyhedron, Facet_id_pmap> Patch_container;
  typedef std::map< Halfedge_handle, Halfedge_handle, Cmp_unik_ad>       Edge_map;
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
  #warning move these to predicates.h
#endif
  template <class Kernel>
  bool  are_triangles_coplanar_same_side(
      const typename Kernel::Point_3& O_prime,const typename Kernel::Point_3& O,
      const typename Kernel::Point_3& P,const typename Kernel::Point_3& Q)
  {
    if ( CGAL::orientation(O_prime, O, P ,Q) != CGAL::COPLANAR ) return false;
    CGAL::Orientation cpl_orient = CGAL::coplanar_orientation(O_prime, O, P, Q);
    CGAL_assertion( cpl_orient != CGAL::COLLINEAR );
    return cpl_orient == CGAL::POSITIVE;
  }


  template <class Nodes_vector, class Vertex_handle>
  bool are_triangles_coplanar_same_side_filtered( int O_prime_index,
                                                 int O_index,
                                                 int P_index,
                                                 int Q_index,
                                                 Vertex_handle P,
                                                 Vertex_handle Q,
                                                 const Nodes_vector& nodes)
  {
    typename Nodes_vector::Protector p;
    try{
      CGAL_USE(p);
      return are_triangles_coplanar_same_side<typename Nodes_vector::Ikernel>(
        nodes.interval_node(O_prime_index),
        nodes.interval_node(O_index),
        P_index == -1 ? nodes.to_interval(get(ppmap,P)): nodes.interval_node(P_index),
        Q_index  == -1 ? nodes.to_interval(get(ppmap,Q)) : nodes.interval_node(Q_index )
      );
    }
    catch(Uncertain_conversion_exception&){
      return are_triangles_coplanar_same_side<typename Nodes_vector::Exact_kernel>(
        nodes.exact_node(O_prime_index),
        nodes.exact_node(O_index),
        P_index == -1 ? nodes.to_exact(get(ppmap,P)): nodes.exact_node(P_index),
        Q_index  == -1 ? nodes.to_exact(get(ppmap,Q)) : nodes.exact_node(Q_index )
      );
    }
  }

  void remove_patches_from_polyhedra(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_to_remove,
    Patch_container& patches_of_P
  )
  {
    for (std::size_t i=patches_to_remove.find_first();
                     i < patches_to_remove.npos; i = patches_to_remove.find_next(i))
    {
      Patch_description& patch=patches_of_P[i];
      internal_IOP::Remove_patch_simplices<Polyhedron>
        modifier( patch.facets,
                  patch.interior_vertices,
                  patch.interior_halfedges,
                  patch.patch_border_halfedges);
      P_ptr->delegate(modifier);
    }
  }

  // function used to remove polylines imported or kept that are incident only
  // to patches not kept for the operation P_ptr is used for storing
  // the result. We look for edges with halfedges both on the border of
  // the mesh. The vertices incident only to such edges should be removed.
  // Here to detect vertices that should be kept, we abuse the fact that
  // the halfedge to be removed and incident to a vertex that should not be
  // removed will still have its next pointer set to a halfedge part of
  // the result.
  void remove_unused_polylines(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_to_remove,
    Patch_container& patches_of_P)
  {
    std::set<Vertex_handle> vertices_to_remove;
    std::set<Halfedge_handle> edges_to_remove;
    for (std::size_t i = patches_to_remove.find_first();
                     i < patches_to_remove.npos;
                     i = patches_to_remove.find_next(i))
    {
      Patch_description& patch=patches_of_P[i];
      BOOST_FOREACH(Halfedge_handle h, patch.patch_border_halfedges)
      {
        if (h->is_border() && h->opposite()->is_border()){
          vertices_to_remove.insert(h->vertex());
          vertices_to_remove.insert(h->opposite()->vertex());
          edges_to_remove.insert( h<h->opposite()?h:h->opposite());
        }
      }
    }

    BOOST_FOREACH(Vertex_handle vh, vertices_to_remove)
    {
      bool to_remove=true;
      BOOST_FOREACH(Halfedge_handle h, halfedges_around_target(vh,*P_ptr))
        if (!h->is_border() || !h->opposite()->is_border())
        {
          to_remove=false;
          // in case the vertex halfedge was one that is going to remove,
          // update it
          set_halfedge(vh, h, *P_ptr);
          break;
        }
      if (to_remove)
        remove_vertex(vh,*P_ptr);
    }
    BOOST_FOREACH(Halfedge_handle hh, edges_to_remove)
      remove_edge(edge(hh,*P_ptr),*P_ptr);
  }

  void disconnect_patches_from_polyhedra(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_to_remove,
    Patch_container& patches_of_P,
    const Edge_map& Phedge_to_Qhedge, //map former patch border halfedge to the equivalent in the other polyhedron
    Edge_map& new_Phedge_to_Qhedge //map new patch border halfedge to the equivalent in the other polyhedron
  )
  {
    for (std::size_t i=patches_to_remove.find_first();
                     i < patches_to_remove.npos; i = patches_to_remove.find_next(i))
    {
      Patch_description& patch=patches_of_P[i];
      std::vector<Halfedge_handle> new_patch_border;
      internal_IOP::Disconnect_patches<Polyhedron>
        modifier( patch.facets,
                  patch.interior_vertices,
                  patch.interior_halfedges,
                  patch.patch_border_halfedges,
                  new_patch_border );
      P_ptr->delegate(modifier);


      CGAL_assertion( new_patch_border.size() ==
                      patch.patch_border_halfedges.size() );

      std::size_t nb_hedges=new_patch_border.size();
      for (std::size_t k=0; k < nb_hedges; ++k){

        CGAL_assertion( patch.patch_border_halfedges[k]->vertex() == new_patch_border[k]->vertex() );
        CGAL_assertion( patch.patch_border_halfedges[k]->opposite()->vertex() == new_patch_border[k]->opposite()->vertex() );
        CGAL_assertion( new_patch_border[k]->is_border_edge() );
        CGAL_assertion( !new_patch_border[k]->is_border() );
        CGAL_assertion( new_patch_border[k]->opposite()->next()->is_border() );
        CGAL_assertion( new_patch_border[k]->opposite()->prev()->is_border() );

        typename Edge_map::const_iterator it_res =
          Phedge_to_Qhedge.find(patch.patch_border_halfedges[k]);
        CGAL_assertion( it_res != Phedge_to_Qhedge.end() );
        CGAL_assertion( it_res->first->vertex()->point() == it_res->second->vertex()->point() );
        CGAL_assertion( it_res->first->opposite()->vertex()->point() == it_res->second->opposite()->vertex()->point() );
        new_Phedge_to_Qhedge[
          patch.patch_border_halfedges[k]==it_res->first
          ? new_patch_border[k]
          : new_patch_border[k]->opposite() ] = it_res->second;
      }

      patch.patch_border_halfedges.swap(new_patch_border);
    }
  }

  template <bool reverse_patch_orientation >
  void append_Q_patches_to_P(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_to_append,
    Patch_container& patches,
    std::map< Halfedge_handle,
            Halfedge_handle,
            internal_IOP::Compare_unik_address<Polyhedron>
          >& Qhedge_to_Phedge
  ){
    #ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning the size of Qhedge_to_Phedge will increase when adding new patches by the size of internal edges. Maybe the use of a copy would be better?
    #endif
    for (std::size_t i=patches_to_append.find_first();
                     i < patches_to_append.npos; i = patches_to_append.find_next(i))
    {
      Patch_description& patch = patches[i];

      internal_IOP::Surface_extension_by_patch_appending<Polyhedron,PolyhedronPointPMap, reverse_patch_orientation>
        modifier(patch.facets, patch.interior_halfedges, patch.patch_border_halfedges, patch.interior_vertices, Qhedge_to_Phedge, ppmap);
      P_ptr->delegate(modifier);
    }

    //CGAL_assertion( Qhedge_to_Phedge.size()==patch_border_halfedges.size() ); // this is not true in case of coplanar patches
    //CGAL_assertion( P_ptr->is_valid() );
  }

  template < class HalfedgeOutputIterator >
  void
  import_polyline(
    Polyhedron& O,
    Halfedge_handle P_first_polyline_hedge,
    Halfedge_handle Q_first_polyline_hedge,
    int nb_segments,
    std::map< Halfedge_handle,
              Halfedge_handle,
              internal_IOP::Compare_unik_address<Polyhedron>
            > & P_to_O_hedge,
    std::map< Halfedge_handle,
              Halfedge_handle,
              internal_IOP::Compare_unik_address<Polyhedron>
            > & Q_to_O_hedge,
    std::map<Vertex_handle, Vertex_handle>& vertex_map,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    HalfedgeOutputIterator output)
  {
    internal_IOP::Import_polyline<Polyhedron, PolyhedronPointPMap, HalfedgeOutputIterator>
      modifier( vertex_map, P_to_O_hedge, Q_to_O_hedge, nb_segments,
                P_first_polyline_hedge, Q_first_polyline_hedge,
                border_halfedges, ppmap, output);

    O.delegate( modifier );
  }

  void compute_border_edge_map(
    const Intersection_polylines& polylines,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    Edge_map& Qhedge_to_Phedge
  ){
    std::size_t nb_polylines = polylines.lengths.size();
    for( std::size_t i=0; i<nb_polylines; ++i)
    {
      if (polylines.to_skip.test(i)) continue;
      Halfedge_handle phedge = polylines.P[i];
      Halfedge_handle qhedge = polylines.Q[i];
      int nb_segments = polylines.lengths[i];

      for (int k=0;k<nb_segments;++k)
      {
        Qhedge_to_Phedge[qhedge]=phedge;
        qhedge=next_marked_halfedge_around_target_vertex(qhedge, patches_of_Q.border_halfedges);
        phedge=next_marked_halfedge_around_target_vertex(phedge, patches_of_P.border_halfedges);
      }
    }
  }

  void compute_inplace_operation(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_of_P_to_keep,
    const boost::dynamic_bitset<>& patches_of_Q_to_import,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    bool reverse_patch_orientation_P,
    bool reverse_patch_orientation_Q,
    Edge_map& Qhedge_to_Phedge
  ){
      //clean up patches not kept
      remove_patches_from_polyhedra(P_ptr, ~patches_of_P_to_keep, patches_of_P);

      if (reverse_patch_orientation_P){
       Polygon_mesh_processing::reverse_face_orientations_of_mesh_with_polylines(*P_ptr);
       // here we need to update the mapping to use the correct border
       // halfedges while appending the patches from Q
       BOOST_FOREACH(typename Edge_map::value_type& v, Qhedge_to_Phedge)
        v.second=v.second->opposite();
     }

      //we import patches from Q
      if ( reverse_patch_orientation_Q )
        append_Q_patches_to_P<true>(P_ptr, patches_of_Q_to_import, patches_of_Q, Qhedge_to_Phedge);
      else
        append_Q_patches_to_P<false>(P_ptr, patches_of_Q_to_import, patches_of_Q, Qhedge_to_Phedge);
  }

  void compute_inplace_operation(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_of_P_to_keep,
    const boost::dynamic_bitset<>& patches_of_Q_to_import,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    bool reverse_patch_orientation_P,
    bool reverse_patch_orientation_Q,
    const Intersection_polylines& polylines
  ){
      Edge_map Qhedge_to_Phedge;
      //maps patch border halfedge from Q to halfedge from P
      compute_border_edge_map(polylines, patches_of_P,
                              patches_of_Q, Qhedge_to_Phedge);

      compute_inplace_operation( P_ptr,
                                 patches_of_P_to_keep, patches_of_Q_to_import,
                                 patches_of_P, patches_of_Q,
                                 reverse_patch_orientation_P, reverse_patch_orientation_Q,
                                 Qhedge_to_Phedge );
    }


#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
  #warning factorize with the code of compute_difference_inplace
#endif
  void compute_inplace_operation_delay_removal_and_insideout(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_of_P_to_keep,
    const boost::dynamic_bitset<>& patches_of_Q_to_import,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    bool reverse_patch_orientation_Q,
    const Intersection_polylines& polylines,
    Edge_map& disconnected_patches_hedge_to_Qhedge
  ){
      Edge_map Qhedge_to_Phedge, Phedge_to_Qhedge;
      //maps patch border halfedge from Q to halfedge from P
      std::size_t nb_polylines = polylines.lengths.size();
      for( std::size_t i=0; i<nb_polylines; ++i)
      {
        Halfedge_handle phedge = polylines.P[i];
        Halfedge_handle qhedge = polylines.Q[i];
        int nb_segments = polylines.lengths[i];

        for (int k=0;;)
        {
          Qhedge_to_Phedge[qhedge]=phedge;
          Phedge_to_Qhedge[phedge]=qhedge;
          if (++k==nb_segments) break;
          qhedge=next_marked_halfedge_around_target_vertex(qhedge, patches_of_Q.border_halfedges);
          phedge=next_marked_halfedge_around_target_vertex(phedge, patches_of_P.border_halfedges);
        }
      }

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
      #warning do not try to disconnect if the patch is isolated? i.e opposite(border_edge_of_patch)->is_border()
#endif
      // disconnect patches inside Q
      // For the patches scheduled to be removed, their patch descriptions
      // in patches_of_P will be updated so that patch_border_halfedges are
      // the newly created halfedges within disconnect_patches_from_polyhedra.
      // Note that disconnected_patches_hedge_to_Qhedge also refers to those halfedges
      //init the map with the previously filled one (needed when reusing patches in two operations)
      disconnected_patches_hedge_to_Qhedge=Phedge_to_Qhedge;
      disconnect_patches_from_polyhedra(P_ptr, ~patches_of_P_to_keep, patches_of_P,
                                        Phedge_to_Qhedge, disconnected_patches_hedge_to_Qhedge);

      //we import patches from Q
      if (reverse_patch_orientation_Q)
        append_Q_patches_to_P<true>(P_ptr, patches_of_Q_to_import, patches_of_Q, Qhedge_to_Phedge);
      else
        append_Q_patches_to_P<false>(P_ptr, patches_of_Q_to_import, patches_of_Q, Qhedge_to_Phedge);
  }

  void remove_disconnected_patches(
    Polyhedron& P,
    Patch_container& patches,
    const boost::dynamic_bitset<>& patches_to_remove)
  {
    for (std::size_t i=patches_to_remove.find_first();
                     i < patches_to_remove.npos;
                     i = patches_to_remove.find_next(i))
    {
      Patch_description& patch=patches[i];
      internal_IOP::Remove_isolated_patch_simplices<Polyhedron>
        modifier( patch.facets, patch.interior_vertices,
          patch.interior_halfedges, patch.patch_border_halfedges
      );
      P.delegate( modifier );
    }
  }

  template < class HalfedgeOutputIterator >
  void fill_new_polyhedron(
    Polyhedron& O, // output
    const boost::dynamic_bitset<>& patches_of_P_to_import,
    const boost::dynamic_bitset<>& patches_of_Q_to_import,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    bool reverse_orientation_of_patches_from_P,
    bool reverse_orientation_of_patches_from_Q,
    const Intersection_polylines& polylines,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    HalfedgeOutputIterator shared_halfedge_output // new shared halfedges
  )
  {
    //add a polyline inside O for each intersection polyline
    std::size_t nb_polylines = polylines.lengths.size();
    std::map<Vertex_handle, Vertex_handle> P_to_O_vertex;
    std::map< Halfedge_handle,
        Halfedge_handle,
        internal_IOP::Compare_unik_address<Polyhedron>
      > P_to_O_hedge, Q_to_O_hedge;
    for (std::size_t i=0; i < nb_polylines; ++i)
      if (!polylines.to_skip.test(i))
        import_polyline(O, polylines.P[i], polylines.Q[i],
                        polylines.lengths[i],
                        P_to_O_hedge, Q_to_O_hedge,
                        P_to_O_vertex, border_halfedges, shared_halfedge_output);

    //import patches of P
    if (reverse_orientation_of_patches_from_P)
      append_Q_patches_to_P<true>(&O, patches_of_P_to_import, patches_of_P, P_to_O_hedge);
    else
      append_Q_patches_to_P<false>(&O, patches_of_P_to_import, patches_of_P, P_to_O_hedge);

    //import patches from Q
    if (reverse_orientation_of_patches_from_Q)
      append_Q_patches_to_P<true>(&O, patches_of_Q_to_import, patches_of_Q, Q_to_O_hedge);
    else
      append_Q_patches_to_P<false>(&O, patches_of_Q_to_import, patches_of_Q, Q_to_O_hedge);
  }

  // detect if a polyline is incident to two patches that won't be imported
  // for the current operation (polylines skipt are always incident to a
  // coplanar patch)
  template <class FacetIdPmap>
  static
  void fill_polylines_to_skip(
    Intersection_polylines& polylines,
    const std::vector<std::size_t>& P_patch_ids,
    const std::vector<std::size_t>& Q_patch_ids,
    const boost::dynamic_bitset<>& patches_of_P_used,
    const boost::dynamic_bitset<>& patches_of_Q_used,
    const FacetIdPmap& P_face_id_pmap,
    const FacetIdPmap& Q_face_id_pmap)
  {
    for (std::size_t i=0;i<polylines.P.size();++i)
    {
      Halfedge_handle h_P = polylines.P[i];
      Halfedge_handle h_Q = polylines.Q[i];
      bool skip_polyline_in_P=true;
      if (!h_P->is_border()){
        std::size_t patch_id = P_patch_ids[ get( P_face_id_pmap, h_P->facet() ) ];
        if (patches_of_P_used.test(patch_id))
          skip_polyline_in_P=false;
      }
      if (skip_polyline_in_P && !h_P->opposite()->is_border()){
        std::size_t patch_id = P_patch_ids[ get( P_face_id_pmap, h_P->opposite()->facet() ) ];
        if (patches_of_P_used.test(patch_id))
          skip_polyline_in_P=false;
      }
      bool skip_polyline_in_Q=true;
      if (!h_Q->is_border()){
        std::size_t patch_id = Q_patch_ids[ get( Q_face_id_pmap, h_Q->facet() ) ];
        if (patches_of_Q_used.test(patch_id))
          skip_polyline_in_Q=false;
      }
      if (skip_polyline_in_Q && !h_Q->opposite()->is_border()){
        std::size_t patch_id = Q_patch_ids[ get( Q_face_id_pmap, h_Q->opposite()->facet() ) ];
        if (patches_of_Q_used.test(patch_id))
          skip_polyline_in_Q=false;
      }

      if (skip_polyline_in_P) polylines.to_skip_in_P.set(i);
      if (skip_polyline_in_Q) polylines.to_skip_in_Q.set(i);
      if (skip_polyline_in_P && skip_polyline_in_Q)
        polylines.to_skip.set(i);
    }
  }

  bool is_dangling_edge(int src_id, int tgt_id,
                        Halfedge_const_handle hedge,
                        const boost::dynamic_bitset<>& is_node_of_degree_one) const
  {
    if ( is_node_of_degree_one.test(src_id) )
    {
      bool res=true;
      Halfedge_const_handle h = hedge->opposite(), start=h;
      do{
        if (h->is_border())
        {
          res = false;
          break;
        }
        h=h->next()->opposite();
      } while(h!=start);
      if (res) return true;
    }
    if ( is_node_of_degree_one.test(tgt_id) )
    {
      Halfedge_const_handle h = hedge, start=h;
      do{
        if (h->is_border())
          return false;
        h=h->next()->opposite();
      } while(h!=start);
      return true;
    }
    return false;
  }

public:

  Polyhedra_output_builder(
    Polyhedron& P,
    Polyhedron& Q,
    cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output_,
    Facet_id_pmap P_facet_id_pmap_,
    Facet_id_pmap Q_facet_id_pmap_,
    PolyhedronPointPMap point_pmap = PolyhedronPointPMap(),
    EdgeMarkPropertyMap edge_pmap = EdgeMarkPropertyMap()
  ) : P_ptr(&P)
    , Q_ptr(&Q)
    , desired_output( desired_output_ )
    , P_facet_id_pmap(P_facet_id_pmap_)
    , Q_facet_id_pmap(Q_facet_id_pmap_)
    , ppmap(point_pmap)
    , edge_mark_pmap(edge_pmap)
    , input_with_coplanar_facets(false)
    , is_P_inside_out( !Polygon_mesh_processing::is_outward_oriented(*P_ptr) )
    , is_Q_inside_out( !Polygon_mesh_processing::is_outward_oriented(*Q_ptr) )
    {}

  Polyhedra_output_builder(
    Polyhedron& P,
    Polyhedron& Q,
    cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output_,
    PolyhedronPointPMap point_pmap = PolyhedronPointPMap(),
    EdgeMarkPropertyMap edge_pmap = EdgeMarkPropertyMap()
  ) : P_ptr(&P)
    , Q_ptr(&Q)
    , desired_output( desired_output_ )
    , ppmap(point_pmap)
    , edge_mark_pmap(edge_pmap)
    , input_with_coplanar_facets(false)
    , is_P_inside_out( !Polygon_mesh_processing::is_outward_oriented(*P_ptr) )
    , is_Q_inside_out( !Polygon_mesh_processing::is_outward_oriented(*Q_ptr) )
    {}

  bool union_valid()        const { return !impossible_operation[P_UNION_Q]; }
  bool intersection_valid() const { return !impossible_operation[P_INTER_Q]; }
  bool P_minus_Q_valid()    const { return !impossible_operation[P_MINUS_Q]; }
  bool Q_minus_P_valid()    const { return !impossible_operation[Q_MINUS_P]; }

  void P_is_inside_out() { is_P_inside_out = true; }
  void Q_is_inside_out() { is_Q_inside_out = true; }

  void input_have_coplanar_facets()
  {
    input_with_coplanar_facets=true;
  }

  template <class Nodes_vector>
  void operator()(
    std::map<Halfedge_const_handle,
             std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    const Nodes_vector& nodes,
    An_edge_per_polyline_map& an_edge_per_polyline,
    const boost::dynamic_bitset<>& is_node_of_degree_one,
    const Poly_to_map_node& /* polyhedron_to_map_node_to_polyhedron_vertex */)
  {
    // first create and fill a map from vertex to node_id
    boost::unordered_map<Vertex_const_handle, int> vertex_to_node_id;
    for(typename std::map<Halfedge_const_handle,
                         std::pair<int,int>,Cmp_unik_ad >::iterator it=border_halfedges.begin(),
                                                                    it_end=border_halfedges.end();
                                                                    it!=it_end;++it)
    {
      vertex_to_node_id[it->first->vertex()]=it->second.second;
      vertex_to_node_id[it->first->opposite()->vertex()]=it->second.first;
    }

    std::size_t num_facets_P = internal::corefinement::init_facet_indices(*P_ptr, P_facet_id_pmap);
    std::size_t num_facets_Q = internal::corefinement::init_facet_indices(*Q_ptr, Q_facet_id_pmap);
    boost::dynamic_bitset<> coplanar_facets_P(num_facets_P, 0);
    boost::dynamic_bitset<> coplanar_facets_Q(num_facets_Q, 0);

    // In the following loop we filter intersection edge that are strictly inside a patch
    // of coplanar facets so that we keep only the edges on the border of the patch.
    // This is not optimal and in an ideal world being able to find the outside edges
    // directly would avoid to compute the intersection of edge/facets inside the patch
    // This loop is done only if the input have some coplanar facets
    typename An_edge_per_polyline_map::iterator epp_it=input_with_coplanar_facets
                                                      ?an_edge_per_polyline.begin()
                                                      :an_edge_per_polyline.end(),
                                                epp_it_end=an_edge_per_polyline.end();
    std::set<Halfedge_const_handle, Cmp_unik_ad> border_edges_to_remove;
    for (;epp_it!=epp_it_end;)
    {
      Halfedge_handle first_hedge  = epp_it->second.first[P_ptr];
      Halfedge_handle first_hedge_opp = first_hedge->opposite();
      Halfedge_handle second_hedge = epp_it->second.first[Q_ptr];
      Halfedge_handle second_hedge_opp = second_hedge->opposite();

      //vertices from P
      // Vertex_handle P1=first_hedge_opp->next()->vertex();
      // Vertex_handle P2=first_hedge->next()->vertex();
      //vertices from Q
      // Vertex_handle Q1=second_hedge_opp->next()->vertex();
      // Vertex_handle Q2=second_hedge->next()->vertex();

      int index_p1=get_node_id(first_hedge_opp->next()->vertex(),vertex_to_node_id);
      int index_p2=get_node_id(first_hedge->next()->vertex(),vertex_to_node_id);
      int index_q1=get_node_id(second_hedge_opp->next()->vertex(),vertex_to_node_id);
      int index_q2=get_node_id(second_hedge->next()->vertex(),vertex_to_node_id);

      // set boolean for the position of P1 wrt to Q1 and Q2
      bool P1_eq_Q1=false, P1_eq_Q2=false;
      if (!first_hedge_opp->is_border() && index_p1!=-1)
      {
        if (!second_hedge_opp->is_border())
          P1_eq_Q1 = index_p1 == index_q1;
        if (!second_hedge->is_border())
          P1_eq_Q2 = index_p1 == index_q2;
      }

      // set boolean for the position of P2 wrt to Q1 and Q2
      bool P2_eq_Q1=false, P2_eq_Q2=false;
      if (!first_hedge->is_border() && index_p2!=-1)
      {
        if (!second_hedge_opp->is_border())
          P2_eq_Q1 = index_p2 == index_q1;
        if (!P1_eq_Q2 && !second_hedge->is_border())
          P2_eq_Q2 = index_p2 == index_q2;
      }

      //mark coplanar facets if any
      if (P1_eq_Q1){
        coplanar_facets_P.set(get(P_facet_id_pmap, first_hedge_opp->facet()));
        coplanar_facets_Q.set(get(Q_facet_id_pmap, second_hedge_opp->facet()));
      }
      if (P1_eq_Q2){
        coplanar_facets_P.set(get(P_facet_id_pmap, first_hedge_opp->facet()));
        coplanar_facets_Q.set(get(Q_facet_id_pmap, second_hedge->facet()));
      }
      if (P2_eq_Q1){
        coplanar_facets_P.set(get(P_facet_id_pmap, first_hedge->facet()));
        coplanar_facets_Q.set(get(Q_facet_id_pmap, second_hedge_opp->facet()));
      }
      if (P2_eq_Q2){
        coplanar_facets_P.set(get(P_facet_id_pmap, first_hedge->facet()));
        coplanar_facets_Q.set(get(Q_facet_id_pmap, second_hedge->facet()));
      }

      // remove the edge if it is in the middle of a coplanar patch
      if ( (P1_eq_Q1 || P1_eq_Q2) && (P2_eq_Q1 || P2_eq_Q2) )
      {
        typename An_edge_per_polyline_map::iterator it_to_rm=epp_it;
        ++epp_it;
        an_edge_per_polyline.erase(it_to_rm);
        border_edges_to_remove.insert(first_hedge);
        border_edges_to_remove.insert(second_hedge);
        #ifdef CGAL_COREFINEMENT_DEBUG
        #warning we need to have the EdgeMarkPropertyMap to unmark intersection hedge
        #endif
      }
      else
        ++epp_it;
    }

    BOOST_FOREACH(Halfedge_const_handle h, border_edges_to_remove)
      border_halfedges.erase(h);

    // (1) Assign a patch id to each facet indicating in which connected
    // component limited by intersection edges of the surface they are.
    internal_IOP::Non_intersection_halfedge<Polyhedron> is_not_marked(border_halfedges);

    CGAL_assertion (P_ptr!=Q_ptr);

    std::vector <std::size_t> P_patch_ids, Q_patch_ids;
    std::vector <std::size_t> P_patch_sizes, Q_patch_sizes;
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning: the extraction of patches could be done by using a seed provided. That is we first make the \
              classification of patches using a representative facet (incident to a marked edge) and once \
              we really need a patch we extract it using this seed. Note that this requires a post-processing \
              if the input polyhedra have several connected component free from intersection \
              this would require to guarantee that there is no connected component not involved in the \
              intersection
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG

    /// \todo I need a property map indicating if an edge is an intersection edge
    std::size_t P_nb_patches = internal::corefinement::
      mark_connected_components_v2(*P_ptr,
                                   is_not_marked,
                                   P_facet_id_pmap,
                                   P_patch_ids,
                                   P_patch_sizes);
    std::size_t Q_nb_patches = internal::corefinement::
      mark_connected_components_v2(*Q_ptr,
                                   is_not_marked,
                                   Q_facet_id_pmap,
                                   Q_patch_ids,
                                   Q_patch_sizes);

    CGAL_assertion ( P_nb_patches==P_patch_sizes.size() );
    CGAL_assertion ( Q_nb_patches==Q_patch_sizes.size() );

    // (2-a) Use the orientation around an edge to classify a patch
    boost::dynamic_bitset<> is_patch_inside_Q(P_nb_patches, false);
    boost::dynamic_bitset<> is_patch_inside_P(Q_nb_patches, false);
    boost::dynamic_bitset<> patch_status_not_set_P(P_nb_patches,true);
    boost::dynamic_bitset<> patch_status_not_set_Q(Q_nb_patches,true);
    boost::dynamic_bitset<> coplanar_patches_of_P(P_nb_patches,false);
    boost::dynamic_bitset<> coplanar_patches_of_Q(Q_nb_patches,false);
    boost::dynamic_bitset<> coplanar_patches_of_P_for_union_and_intersection(P_nb_patches,false);
    boost::dynamic_bitset<> coplanar_patches_of_Q_for_union_and_intersection(Q_nb_patches,false);

    for (typename An_edge_per_polyline_map::const_iterator it=an_edge_per_polyline.begin();it!=an_edge_per_polyline.end();++it)
    {
      CGAL_assertion(it->second.first.size()==2);
      //orientation of faces around the edge (to be sure we can do it)
      std::pair<int,int> indices = it->first;
      //const std::pair<bool,int>& polyline_info=it->second.second;

      //get the two halfedges incident to the edge [indices.first,indices.second]
      Halfedge_handle first_hedge  = it->second.first.find(P_ptr)->second;
      Halfedge_handle second_hedge = it->second.first.find(Q_ptr)->second;

      CGAL_assertion(nodes[indices.second]==get(ppmap,first_hedge->vertex()));
      CGAL_assertion(nodes[indices.first]==get(ppmap,first_hedge->opposite()->vertex()));
      CGAL_assertion(nodes[indices.second]==get(ppmap,second_hedge->vertex()));
      CGAL_assertion(nodes[indices.first]==get(ppmap,second_hedge->opposite()->vertex()));

      //different handling depending on the number of incident triangles to the edge.
      //After sewing there are two,three or four volumes if there are two,three or four incident triangles respectively
      if ( first_hedge->is_border_edge() ){
        if (second_hedge->is_border_edge() )
        {
          if ( first_hedge->is_border() != second_hedge->is_border() )
          {
            //No restriction at this level
            std::size_t patch_id_P =
              P_patch_ids[ get( P_facet_id_pmap, first_hedge->is_border()
                                                 ? first_hedge->opposite()->facet()
                                                 : first_hedge->facet() ) ];
            std::size_t patch_id_Q =
              Q_patch_ids[ get( Q_facet_id_pmap, second_hedge->is_border()
                                                 ? second_hedge->opposite()->facet()
                                                 : second_hedge->facet() ) ];
            patch_status_not_set_P.reset(patch_id_P);
            patch_status_not_set_Q.reset(patch_id_Q);
          }
          else
          {
            //Nothing allowed
            impossible_operation.set();
            return;
          }
        }
        else
        {
          //Ambiguous, we can do nothing
          impossible_operation.set();
          return;
        }
      }
      else
        if ( second_hedge->is_border_edge() )
        {
          //Ambiguous, we do nothing
          impossible_operation.set();
          return;
        }
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
          int index_p1=get_node_id(first_hedge->opposite()->next()->vertex(),vertex_to_node_id);
          int index_p2=get_node_id(first_hedge->next()->vertex(),vertex_to_node_id);
          int index_q1=get_node_id(second_hedge->opposite()->next()->vertex(),vertex_to_node_id);
          int index_q2=get_node_id(second_hedge->next()->vertex(),vertex_to_node_id);

          std::size_t patch_id_p1=P_patch_ids[ get(P_facet_id_pmap, first_hedge->opposite()->facet()) ];
          std::size_t patch_id_p2=P_patch_ids[ get(P_facet_id_pmap, first_hedge->facet()) ];
          std::size_t patch_id_q1=Q_patch_ids[ get(Q_facet_id_pmap, second_hedge->opposite()->facet()) ];
          std::size_t patch_id_q2=Q_patch_ids[ get(Q_facet_id_pmap, second_hedge->facet()) ];

          //indicates that patch status will be updated
          patch_status_not_set_P.reset(patch_id_p1);
          patch_status_not_set_P.reset(patch_id_p2);
          patch_status_not_set_Q.reset(patch_id_q1);
          patch_status_not_set_Q.reset(patch_id_q2);

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
          #warning: Factorize the orientation predicates.
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
          // handle case of coplanar facets
          // We choose that a coplanar patch is classified like the other incident patch since they bound the same volume.
          if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p1,index_q1,P1,Q1,nodes) ) //P1==Q1
          {
            coplanar_patches_of_P.set(patch_id_p1);
            coplanar_patches_of_Q.set(patch_id_q1);
            coplanar_patches_of_P_for_union_and_intersection.set(patch_id_p1);
            coplanar_patches_of_Q_for_union_and_intersection.set(patch_id_q1);

            CGAL_assertion( !are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p2,index_q2,P2,Q2,nodes) );

            bool Q2_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);
            if ( Q2_is_between_P1P2 ) is_patch_inside_P.set(patch_id_q2); //case 1
            else is_patch_inside_Q.set(patch_id_p2); //case 2
            continue;
          }
          else{
            if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p1,index_q2,P1,Q2,nodes) ) //P1==Q2
            {
              CGAL_assertion( index_p1!=index_p2 || index_p1==-1 );
              coplanar_patches_of_P.set(patch_id_p1);
              coplanar_patches_of_Q.set(patch_id_q2);
              bool Q1_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
              if ( Q1_is_between_P1P2 )
              { // case 3
                is_patch_inside_P.set(patch_id_q1);
                is_patch_inside_Q.set(patch_id_p2);
              } //else case 4
              continue;
            }
            else
            {
              if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p2,index_q1,P2,Q1,nodes) ) //P2==Q1
              {
                coplanar_patches_of_P.set(patch_id_p2);
                coplanar_patches_of_Q.set(patch_id_q1);
                bool Q2_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);
                if ( Q2_is_between_P1P2 )
                {  //case 5
                  is_patch_inside_P.set(patch_id_q2);
                  is_patch_inside_Q.set(patch_id_p1);
                } // else case 6
                continue;
              }
              else{
                if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p2,index_q2,P2,Q2,nodes) ) //P2==Q2
                {
                  coplanar_patches_of_P.set(patch_id_p2);
                  coplanar_patches_of_Q.set(patch_id_q2);
                  coplanar_patches_of_P_for_union_and_intersection.set(patch_id_p2);
                  coplanar_patches_of_Q_for_union_and_intersection.set(patch_id_q2);
                  bool Q1_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
                  if ( Q1_is_between_P1P2 ) is_patch_inside_P.set(patch_id_q1);  //case 7
                  else is_patch_inside_Q.set(patch_id_p1); //case 8
                  continue;
                }
              }
            }
          }
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
          #warning At some point we should have a check if a patch status is already set, what we do is consistant otherwise --> ambiguous
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG

          CGAL_assertion(
              ( index_p1 == -1 ? nodes.to_exact(get(ppmap,P1)): nodes.exact_node(index_p1) ) !=
              ( index_q1 == -1 ? nodes.to_exact(get(ppmap,Q1)): nodes.exact_node(index_q1) )
          &&
              ( index_p2 == -1 ? nodes.to_exact(get(ppmap,P2)): nodes.exact_node(index_p2) ) !=
              ( index_q1 == -1 ? nodes.to_exact(get(ppmap,Q1)): nodes.exact_node(index_q1) )
          &&
              ( index_p1 == -1 ? nodes.to_exact(get(ppmap,P1)): nodes.exact_node(index_p1) ) !=
              ( index_q2 == -1 ? nodes.to_exact(get(ppmap,Q2)): nodes.exact_node(index_q2) )
          &&
              ( index_p2 == -1 ? nodes.to_exact(get(ppmap,P2)): nodes.exact_node(index_p2) ) !=
              ( index_q2 == -1 ? nodes.to_exact(get(ppmap,Q2)): nodes.exact_node(index_q2) )
          );

          bool Q1_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
          bool Q2_is_between_P1P2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);

          if ( Q1_is_between_P1P2 ){
            is_patch_inside_P.set(patch_id_q1);
            if( Q2_is_between_P1P2 )
            {
              is_patch_inside_P.set(patch_id_q2);
              bool P1_is_between_Q1Q2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes,ppmap);
              if (!P1_is_between_Q1Q2){
                // case (a4)
                // poly_first  - poly_second            = P1Q1 U Q2P2
                // poly_second - poly_first             = {0}
                // poly_first \cap poly_second          = Q1Q2
                // opposite( poly_first U poly_second ) = P2P1
                impossible_operation.set(P_MINUS_Q); // P-Q is non-manifold
              }
              else{
                // case (b4)
                // poly_first  - poly_second            = Q2Q1
                // poly_second - poly_first             = P2P1
                // poly_first \cap poly_second          = P1Q2 U Q1P2
                // opposite( poly_first U poly_second ) = {O}
                is_patch_inside_Q.set(patch_id_p1);
                is_patch_inside_Q.set(patch_id_p2);
                impossible_operation.set(P_INTER_Q); // P \cap Q is non-manifold
              }
            }
            else
            {
              //case (c4)
              // poly_first  - poly_second            = P1Q1
              // poly_second - poly_first             = P2Q2
              // poly_first \cap poly_second          = Q1P2
              // opposite( poly_first U poly_second ) = Q2P1
              if ( is_dangling_edge(indices.first, indices.second, first_hedge, is_node_of_degree_one) ||
                   is_dangling_edge(indices.first, indices.second, second_hedge, is_node_of_degree_one) )
              {
                impossible_operation.set();
                return;
              }
              is_patch_inside_Q.set(patch_id_p2);
            }
          }
          else
          {
            if( Q2_is_between_P1P2 )
            {
              //case (d4)
              // poly_first  - poly_second            = Q2P2
              // poly_second - poly_first             = Q1P1
              // poly_first \cap poly_second          = P1Q2
              // opposite( poly_first U poly_second ) = P2Q1
              if ( is_dangling_edge(indices.first, indices.second, first_hedge, is_node_of_degree_one) ||
                   is_dangling_edge(indices.first, indices.second, second_hedge, is_node_of_degree_one) )
              {
                impossible_operation.set();
                return;
              }
              is_patch_inside_P.set(patch_id_q2);
              is_patch_inside_Q.set(patch_id_p1);
            }
            else
            {
              bool P1_is_between_Q1Q2 = OOP::sorted_around_edge_filtered(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes,ppmap);
              if (!P1_is_between_Q1Q2){
                //case (e4)
                // poly_first  - poly_second            = P1P2
                // poly_second - poly_first             = Q1Q2
                // poly_first \cap poly_second          = {0}
                // opposite( poly_first U poly_second ) = P2Q1 U Q2P1
                impossible_operation.set(P_UNION_Q); // P U Q is non-manifold
              }
              else{
                //case (f4)
                is_patch_inside_Q.set(patch_id_p1);
                is_patch_inside_Q.set(patch_id_p2);
                // poly_first  - poly_second            = {0}
                // poly_second - poly_first             = Q1P1 U P2Q2
                // poly_first \cap poly_second          = P1P2
                // opposite( poly_first U poly_second ) = Q2Q1
                impossible_operation.set(Q_MINUS_P); // Q - P is non-manifold
              }
            }
          }
        }
    }


///////////////////////////////////////////////
//////////////////////END//////////////////////
///////////////////////////////////////////////

    // (2-b) Classify isolated surface patches wrt the other support polyhedron
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning this does not work with open polyhedra
    #warning this should not be done if we have surfaces with boundaries!!! ask the user a flag?
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    typedef Side_of_triangle_mesh<Polyhedron, Kernel, PolyhedronPointPMap> Inside_poly_test;

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning stop using next_marked_halfedge_around_target_vertex and create lists of halfedges instead?
#endif

    if ( patch_status_not_set_P.any() )
    {
      CGAL::Bounded_side in_Q = is_Q_inside_out ? ON_UNBOUNDED_SIDE : ON_BOUNDED_SIDE;

      Inside_poly_test inside_Q(*Q_ptr, ppmap);
      for (typename Polyhedron::Face_iterator fit=P_ptr->facets_begin(),
                                              fit_end=P_ptr->facets_end();
                                              fit!=fit_end; ++fit)
      {
        std::size_t patch_id=P_patch_ids[ get(P_facet_id_pmap, fit) ];
        if ( patch_status_not_set_P.test( patch_id ) )
        {
          patch_status_not_set_P.reset( patch_id );
          CGAL::Bounded_side position = inside_Q( get(ppmap,fit->halfedge()->vertex()) );
          if ( position == in_Q )
            is_patch_inside_Q.set(patch_id);
          else
            if ( position == ON_BOUNDARY )
            {
              if (coplanar_facets_P.test(get(P_facet_id_pmap, fit)))
              {
                coplanar_patches_of_P.set(patch_id);
                coplanar_patches_of_P_for_union_and_intersection.set(patch_id);
              }
              else
              {
                Vertex_handle vn = fit->halfedge()->opposite()->vertex();
                Bounded_side other_position = inside_Q( get(ppmap, vn) );
                if (other_position==ON_BOUNDARY)
                {
                  // \todo improve this part which is not robust with a kernel with inexact constructions.
                  other_position = inside_Q(midpoint(get(ppmap, vn),
                                                     get(ppmap, fit->halfedge()->vertex()) ));
                }
                if ( other_position == in_Q )
                  is_patch_inside_Q.set(patch_id);
              }
            }

          if ( patch_status_not_set_P.none() ) break;
        }
      }
    }

    if ( patch_status_not_set_Q.any() )
    {
      CGAL::Bounded_side in_P = is_P_inside_out ? ON_UNBOUNDED_SIDE : ON_BOUNDED_SIDE;

      Inside_poly_test inside_P(*P_ptr, ppmap);
      for (typename Polyhedron::Face_iterator fit=Q_ptr->facets_begin(),
                                              fit_end=Q_ptr->facets_end();
                                              fit!=fit_end; ++fit)
      {
        std::size_t patch_id=Q_patch_ids[ get(Q_facet_id_pmap, fit) ];
        if ( patch_status_not_set_Q.test( patch_id ) )
        {
          patch_status_not_set_Q.reset( patch_id );
          Bounded_side position = inside_P( get(ppmap,fit->halfedge()->vertex()) );
          if ( position == in_P )
            is_patch_inside_P.set(patch_id);
          else
            if ( position == ON_BOUNDARY )
            {
              if (coplanar_facets_Q.test(get(Q_facet_id_pmap, fit)))
              {
                coplanar_patches_of_Q.set(patch_id);
                coplanar_patches_of_Q_for_union_and_intersection.set(patch_id);
              }
              else
              {
                Vertex_handle vn = fit->halfedge()->opposite()->vertex();
                Bounded_side other_position = inside_P( get(ppmap, vn) );
                if (other_position==ON_BOUNDARY)
                {
                  // \todo improve this part which is not robust with a kernel with inexact constructions.
                  other_position = inside_P(midpoint(get(ppmap, vn),
                                                     get(ppmap, fit->halfedge()->vertex()) ));
                }
                if ( other_position == in_P )
                  is_patch_inside_P.set(patch_id);
              }
            }
          if ( patch_status_not_set_Q.none() ) break;
        }
      }
    }

    /// \todo There might be some patches for which the status is unknown, we need to use
    ///       for example the centroid of the facet to decide, but this would require an exact
    ///       inside-polyhedron object or to work with points at endpoints of the intervals.
    ///       See SL_TMP_ASSERT
    //to maintain a polyhedron halfedge on each polyline + pair<bool,int>
    //with first = "is the key (pair<int,int>) was reversed?" and second is the number of edges +1 in the polyline
    //typedef std::map< std::pair<int,int>, std::pair< std::map<Polyhedron*,Halfedge_handle>,std::pair<bool,int> > > An_edge_per_polyline_map;

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning add a mechanism to handle the patches independantly (for example calculating the volume without building the polyhedron) \
             This can be done by using a functor to which we give the bitset, the polyhedra and one facet per patch?
#endif // CGAL_COREFINEMENT_POLYHEDRA_DEBUG

    #ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "is_patch_inside_Q " <<  is_patch_inside_Q << "\n";
    std::cout << "is_patch_inside_P " << is_patch_inside_P << "\n";
    std::cout << "coplanar_patches_of_P " << coplanar_patches_of_P << "\n";
    std::cout << "coplanar_patches_of_Q " << coplanar_patches_of_Q << "\n";
    std::cout << "coplanar_patches_of_P_for_union_and_intersection " << coplanar_patches_of_P_for_union_and_intersection << "\n";
    std::cout << "coplanar_patches_of_Q_for_union_and_intersection " << coplanar_patches_of_Q_for_union_and_intersection << "\n";
    std::cout << "Size of patches of P: ";
    std::copy(P_patch_sizes.rbegin(), P_patch_sizes.rend(), std::ostream_iterator<std::size_t>(std::cout," ") );
    std::cout << "\n";
    std::cout << "Size of patches of Q: ";
    std::copy(Q_patch_sizes.rbegin(), Q_patch_sizes.rend(), std::ostream_iterator<std::size_t>(std::cout," ") );
    std::cout << "\n";
    #endif

    //backup an halfedge per polyline
    std::vector <Halfedge_handle> P_polylines, Q_polylines;
    std::vector<int> polyline_lengths;

    for (typename An_edge_per_polyline_map::const_iterator
          it=an_edge_per_polyline.begin(),
          it_end=an_edge_per_polyline.end();
          it!=it_end;++it)
    {
      const std::pair<bool,int>& polyline_info=it->second.second;

      Halfedge_handle qhedge = it->second.first.find( Q_ptr )->second;
      Halfedge_handle phedge = it->second.first.find( P_ptr )->second;

      if( polyline_info.first ){
        phedge=phedge->opposite();
        qhedge=qhedge->opposite();
      }

      P_polylines.push_back(phedge);
      Q_polylines.push_back(qhedge);
      polyline_lengths.push_back(polyline_info.second+1);
    }

    //store the patch description in a container to avoid recomputing it several times
    Patch_container patches_of_P( P_ptr, P_patch_ids, P_facet_id_pmap, border_halfedges, P_nb_patches),
                    patches_of_Q( Q_ptr, Q_patch_ids, Q_facet_id_pmap, border_halfedges, Q_nb_patches);

    // for each boolean operation, define two bitsets of patches contributing
    // to the result
    std::vector< boost::dynamic_bitset<> > patches_of_P_used(4);
    std::vector< boost::dynamic_bitset<> > patches_of_Q_used(4);

    /// handle the bitset for the union
    if ( !impossible_operation.test(P_UNION_Q) && desired_output[P_UNION_Q] )
    {
      //define patches to import from P
      patches_of_P_used[P_UNION_Q] = ~is_patch_inside_Q - coplanar_patches_of_P;
      //define patches to import from Q
      patches_of_Q_used[P_UNION_Q] = ~is_patch_inside_P - coplanar_patches_of_Q;
      //handle coplanar patches
      if (coplanar_patches_of_P.any())
      {
        if (desired_output[P_UNION_Q]==Q_ptr)
          patches_of_Q_used[P_UNION_Q] |= coplanar_patches_of_Q_for_union_and_intersection;
        else
          patches_of_P_used[P_UNION_Q] |= coplanar_patches_of_P_for_union_and_intersection;
      }
    }

    /// handle the bitset for the intersection
    if ( !impossible_operation.test(P_INTER_Q) && desired_output[P_INTER_Q] )
    {
      //define patches to import from P
      patches_of_P_used[P_INTER_Q] = is_patch_inside_Q;
      //define patches to import from Q
      patches_of_Q_used[P_INTER_Q] = is_patch_inside_P;
      //handle coplanar patches
      if (coplanar_patches_of_P.any())
      {
        if (desired_output[P_INTER_Q]==Q_ptr)
          patches_of_Q_used[P_INTER_Q] |= coplanar_patches_of_Q_for_union_and_intersection;
        else
          patches_of_P_used[P_INTER_Q] |= coplanar_patches_of_P_for_union_and_intersection;
      }
    }

    /// handle the bitset for P-Q
    if ( !impossible_operation.test(P_MINUS_Q) && desired_output[P_MINUS_Q] )
    {
      //define patches to import from P
      patches_of_P_used[P_MINUS_Q] = (~is_patch_inside_Q - coplanar_patches_of_P);
      //define patches to import from Q
      patches_of_Q_used[P_MINUS_Q] = is_patch_inside_P;
      //handle coplanar patches
      if (coplanar_patches_of_P.any())
      {
        if (desired_output[P_MINUS_Q]==Q_ptr)
          patches_of_Q_used[P_MINUS_Q] |= ~coplanar_patches_of_Q_for_union_and_intersection & coplanar_patches_of_Q;
        else
          patches_of_P_used[P_MINUS_Q] |= ~coplanar_patches_of_P_for_union_and_intersection & coplanar_patches_of_P;
      }
    }

    /// handle the bitset for Q-P
    if ( !impossible_operation.test(Q_MINUS_P) && desired_output[Q_MINUS_P] )
    {
      //define patches to import from P
      patches_of_P_used[Q_MINUS_P] = is_patch_inside_Q;
      //define patches to import from Q
      patches_of_Q_used[Q_MINUS_P] = ~is_patch_inside_P - coplanar_patches_of_Q;
      //handle coplanar patches
      if (coplanar_patches_of_P.any())
      {
        if (desired_output[Q_MINUS_P]==Q_ptr)
          patches_of_Q_used[Q_MINUS_P] |= ~coplanar_patches_of_Q_for_union_and_intersection & coplanar_patches_of_Q;
        else
          patches_of_P_used[Q_MINUS_P] |= ~coplanar_patches_of_P_for_union_and_intersection & coplanar_patches_of_P;
      }
    }

    #ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "patches_of_P_used[P_UNION_Q] " << patches_of_P_used[P_UNION_Q] << "\n";
    std::cout << "patches_of_Q_used[P_UNION_Q] " << patches_of_Q_used[P_UNION_Q] << "\n";
    std::cout << "patches_of_P_used[P_INTER_Q] " << patches_of_P_used[P_INTER_Q] << "\n";
    std::cout << "patches_of_Q_used[P_INTER_Q] " << patches_of_Q_used[P_INTER_Q] << "\n";
    std::cout << "patches_of_P_used[P_MINUS_Q] " << patches_of_P_used[P_MINUS_Q] << "\n";
    std::cout << "patches_of_Q_used[P_MINUS_Q] " << patches_of_Q_used[P_MINUS_Q] << "\n";
    std::cout << "patches_of_P_used[Q_MINUS_P] " << patches_of_P_used[Q_MINUS_P] << "\n";
    std::cout << "patches_of_Q_used[Q_MINUS_P] " << patches_of_Q_used[Q_MINUS_P] << "\n";
    #endif // CGAL_COREFINEMENT_DEBUG
    // Schedule the order in which the different boolean operations should be
    // done. First operations are those filling polyhedra different
    // from P and Q, then the one modifying P and finally the one
    // modifying Q.
    std::vector<Boolean_operation> out_of_place_operations;
    Boolean_operation inplace_operation_P=NONE, inplace_operation_Q=NONE;
    for (int i=0;i<4;++i)
    {
      Boolean_operation operation=enum_cast<Boolean_operation>(i);

      if (!desired_output[operation] || impossible_operation.test(operation))
        continue;

      if (desired_output[operation]==P_ptr)
        inplace_operation_P=operation;
      else
        if (desired_output[operation]==Q_ptr)
          inplace_operation_Q=operation;
        else
          out_of_place_operations.push_back(operation);
    }

    /// first handle operations in a polyhedron that is neither P nor Q
    BOOST_FOREACH(Boolean_operation operation, out_of_place_operations)
    {
      Polyhedron* ouput_ptr = *desired_output[operation];
      CGAL_assertion(P_ptr!=ouput_ptr && Q_ptr!=ouput_ptr);

      Intersection_polylines polylines(P_polylines, Q_polylines, polyline_lengths);
      // skip the import of polylines only incident to patch(es)
      // not used by the current operation
      fill_polylines_to_skip(
        polylines, P_patch_ids, Q_patch_ids,
        patches_of_P_used[operation], patches_of_Q_used[operation],
        P_facet_id_pmap, Q_facet_id_pmap
      );

      std::vector<Halfedge_handle> shared_halfedges;
      fill_new_polyhedron(
        *ouput_ptr,
        patches_of_P_used[operation], patches_of_Q_used[operation],
        patches_of_P, patches_of_Q,
        operation == Q_MINUS_P, operation == P_MINUS_Q,
        polylines,
        border_halfedges,
        std::back_inserter(shared_halfedges)
      );
      BOOST_FOREACH(Halfedge_handle h, shared_halfedges){
        put(edge_mark_pmap, std::make_pair(h,ouput_ptr),true);
        put(edge_mark_pmap, std::make_pair(h->opposite(),ouput_ptr),true);
      }
    }

    Edge_map disconnected_patches_hedge_to_Qhedge;

    /// handle the operations updating P and/or Q
    if ( inplace_operation_P!=NONE )
    {
      CGAL_assertion( *desired_output[inplace_operation_P] == P_ptr );

      if ( inplace_operation_Q!=NONE)
      {
        // operation in P with removal (and optinally inside-out) delayed
        // First backup the border edges of patches to be used
        Patch_container tmp_patches_of_P(P_ptr,
          patches_of_P.patch_ids,
          patches_of_P.facet_id_pmap,
          patches_of_P.border_halfedges,
          patches_of_P.patches.size());
        boost::dynamic_bitset<> patches_of_P_removed = ~patches_of_P_used[inplace_operation_P];
        for (std::size_t i = patches_of_P_removed.find_first();
                         i < patches_of_P_removed.npos;
                         i = patches_of_P_removed.find_next(i))
        {
          // we are only interested by patch border halfedges so
          // squeeze the auto-filling mechanism
          tmp_patches_of_P.patches[i].is_initialized=true;
          tmp_patches_of_P.patches[i].patch_border_halfedges=
            patches_of_P[i].patch_border_halfedges;
        }

        Intersection_polylines polylines_in_P(
          P_polylines, Q_polylines, polyline_lengths);
        Intersection_polylines polylines_in_Q=polylines_in_P;
        fill_polylines_to_skip(
          polylines_in_P, P_patch_ids, Q_patch_ids,
          patches_of_P_used[inplace_operation_P],
          patches_of_Q_used[inplace_operation_P],
          P_facet_id_pmap, Q_facet_id_pmap);
        fill_polylines_to_skip(
          polylines_in_Q, P_patch_ids, Q_patch_ids,
          patches_of_P_used[inplace_operation_Q],
          patches_of_Q_used[inplace_operation_Q],
          P_facet_id_pmap, Q_facet_id_pmap);
        // force the initialization of the patches of P used
        // for the operation in Q before P is modified
        for (std::size_t i=patches_of_P_used[inplace_operation_Q].find_first();
                         i < patches_of_P_used[inplace_operation_Q].npos;
                         i = patches_of_P_used[inplace_operation_Q].find_next(i))
        {
          patches_of_P[i];
        }
        // Operation in P: disconnect patches not use and append the one from Q
        compute_inplace_operation_delay_removal_and_insideout(
          P_ptr,
          patches_of_P_used[inplace_operation_P], patches_of_Q_used[inplace_operation_P],
          patches_of_P, patches_of_Q,
          inplace_operation_P == P_MINUS_Q || inplace_operation_P == Q_MINUS_P,
          polylines_in_P, disconnected_patches_hedge_to_Qhedge);
        // Operation in Q: discard patches and append the one from Q
        CGAL_assertion( *desired_output[inplace_operation_Q] == Q_ptr );
        compute_inplace_operation( Q_ptr,
                                   patches_of_Q_used[inplace_operation_Q],
                                   patches_of_P_used[inplace_operation_Q],
                                   patches_of_Q, patches_of_P,
                                   inplace_operation_Q==P_MINUS_Q,
                                   inplace_operation_Q==Q_MINUS_P,
                                   disconnected_patches_hedge_to_Qhedge);
        // remove polylines only on the border of patches not kept in Q
        if (polylines_in_Q.to_skip.any())
          remove_unused_polylines(Q_ptr,
                                  ~patches_of_Q_used[inplace_operation_Q],
                                  patches_of_Q);
        // now remove patches temporarily kept in P
         remove_disconnected_patches(*P_ptr, patches_of_P, patches_of_P_removed);
        // remove polylines only on the border of patches not kept in P
        if (polylines_in_P.to_skip.any())
          remove_unused_polylines(P_ptr,
                                  ~patches_of_P_used[inplace_operation_P],
                                  tmp_patches_of_P);
         // finally reverse orientation of P if needed
         if (inplace_operation_P == Q_MINUS_P)
           CGAL::Polygon_mesh_processing::reverse_face_orientations(*P_ptr);
      }
      else{
        /// handle the operation updating only P
        CGAL_assertion( *desired_output[inplace_operation_P] == P_ptr );
        Intersection_polylines polylines(
          P_polylines, Q_polylines, polyline_lengths);
        fill_polylines_to_skip(
          polylines, P_patch_ids, Q_patch_ids,
          patches_of_P_used[inplace_operation_P],
          patches_of_Q_used[inplace_operation_P],
          P_facet_id_pmap, Q_facet_id_pmap
        );

        compute_inplace_operation(
          P_ptr,
          patches_of_P_used[inplace_operation_P],
          patches_of_Q_used[inplace_operation_P],
          patches_of_P, patches_of_Q,
          inplace_operation_P == Q_MINUS_P,
          inplace_operation_P == P_MINUS_Q,
          polylines
        );
        // remove polylines only on the border of patches not kept
        if (polylines.to_skip.any())
          remove_unused_polylines(P_ptr,
                                  ~patches_of_P_used[inplace_operation_P],
                                  patches_of_P);
      }
    }
    else
      if ( inplace_operation_Q!=NONE )
      {
        /// handle the operation updating only Q
        CGAL_assertion( *desired_output[inplace_operation_Q] == Q_ptr );
        Intersection_polylines polylines(
          Q_polylines, P_polylines, polyline_lengths);
        fill_polylines_to_skip(
          polylines, Q_patch_ids, P_patch_ids,
          patches_of_Q_used[inplace_operation_Q],
          patches_of_P_used[inplace_operation_Q],
          Q_facet_id_pmap, P_facet_id_pmap
        );

        compute_inplace_operation( Q_ptr,
                                   patches_of_Q_used[inplace_operation_Q],
                                   patches_of_P_used[inplace_operation_Q],
                                   patches_of_Q, patches_of_P,
                                   inplace_operation_Q==P_MINUS_Q,
                                   inplace_operation_Q==Q_MINUS_P,
                                   polylines);

        // remove polylines only on the border of patches not kept
        if (polylines.to_skip.any())
          remove_unused_polylines(Q_ptr,
                                  ~patches_of_Q_used[inplace_operation_Q],
                                  patches_of_Q);
      }
  }
};

} } // end of namespace CGAL::Corefinement

#include <CGAL/enable_warnings.h>

#endif // CGAL_INTERNAL_COREFINEMENT_POLYHEDRA_OUTPUT_BUILDER_H
