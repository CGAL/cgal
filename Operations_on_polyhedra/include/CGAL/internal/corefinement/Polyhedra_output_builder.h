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

#ifndef CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H
#define CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H

#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>
#include <CGAL/Point_inside_polyhedron_3.h>
#include <CGAL/Triangle_accessor_with_ppmap_3.h>
#include <CGAL/internal/corefinement/predicates.h>
#include <CGAL/internal/corefinement/utils.h>
#include <CGAL/array.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>

namespace CGAL{

namespace internal_IOP{


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
  while(is_marked.find(next)==is_marked.end());
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

template < class Polyhedron, class PolyhedronPointPMap >
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

  Import_polyline(
    Vertex_map& vertex_map_,
    Hedge_map&  P_to_O_hedge_,
    Hedge_map&  Q_to_O_hedge_,
    int nb_segments_,
    Halfedge_handle P_first_halfedge_,
    Halfedge_handle Q_first_halfedge_,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges_,
    PolyhedronPointPMap ppmap_)
  : vertex_map( vertex_map_ )
  , P_to_O_hedge( P_to_O_hedge_ )
  , Q_to_O_hedge( Q_to_O_hedge_ )
  , nb_segments( nb_segments_ )
  , P_first_halfedge( P_first_halfedge_ )
  , Q_first_halfedge( Q_first_halfedge_ )
  , border_halfedges( border_halfedges_ )
  , ppmap( ppmap_ )
  {}

  void operator()(HDS& hds)
  {
    HalfedgeDS_decorator<HDS> decorator(hds);

    typename HDS::Halfedge dummy_hedge; //dummy default constructed halfedge
    O_first_halfedge=hds.edges_push_back(dummy_hedge, dummy_hedge);

    //make sure the first vertex does not already exist
    Vertex_handle source;
    std::pair< typename std::map< Vertex_handle, Vertex_handle >::iterator, bool > insert_res=
      vertex_map.insert( std::make_pair( P_first_halfedge->opposite()->vertex(), source ) );


    if( insert_res.second )
    {
      source = decorator.vertices_push_back( *(P_first_halfedge->opposite()->vertex()) );
      put(ppmap, source, get(ppmap, P_first_halfedge->opposite()->vertex() ) );
      decorator.set_vertex_halfedge(source, O_first_halfedge->opposite());
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
        put(ppmap, target, get(ppmap, P_first_halfedge->vertex() ) );
        decorator.set_vertex_halfedge(target, O_first_halfedge);
        insert_res.first->second = target;
      }
      else
        target = insert_res.first->second;
    }
    else{
      target=decorator.vertices_push_back( *(P_first_halfedge->vertex()) );
      put(ppmap, target, get(ppmap, P_first_halfedge->vertex() ) );
      decorator.set_vertex_halfedge(target, O_first_halfedge);
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
      //get the new edge
      Halfedge_handle P_hedge = next_marked_halfedge_around_target_vertex(P_previous, border_halfedges);
      Halfedge_handle Q_hedge = next_marked_halfedge_around_target_vertex(Q_previous, border_halfedges);
      //if this is the final segment, only create a target vertex if it does not exist
      if (i!=nb_segments-1)
      {
        target=decorator.vertices_push_back( *(P_hedge->vertex()) );
        put(ppmap, target, get(ppmap, P_hedge->vertex() ) );
        decorator.set_vertex_halfedge(target, O_hedge);
      }
      else{
        std::pair< typename std::map< Vertex_handle, Vertex_handle >::iterator, bool > insert_res=
          vertex_map.insert( std::make_pair(P_hedge->vertex(), Vertex_handle()) );
        if (insert_res.second)
        {
          target=decorator.vertices_push_back( *(P_hedge->vertex()) );
          put(ppmap, target, get(ppmap, P_hedge->vertex() ) );
          decorator.set_vertex_halfedge(target, O_hedge);
          insert_res.first->second = target;
        }
        else
          target = insert_res.first->second;
      }

      //set next/prev relationships
      O_previous->HBase::set_next( O_hedge );
      O_hedge->opposite()->HBase::set_next( O_previous->opposite() );
      decorator.set_prev(O_hedge, O_previous);
      decorator.set_prev(O_previous->opposite(), O_hedge->opposite());

      decorator.set_vertex(O_hedge, target);
      decorator.set_vertex(O_hedge->opposite(), source);

      O_previous=O_hedge;
      P_previous = P_hedge;
      Q_previous = Q_hedge;
      source = target;

      P_to_O_hedge.insert( std::make_pair(P_previous, O_previous) );
      Q_to_O_hedge.insert( std::make_pair(Q_previous, O_previous) );
    }

    //update next/prev relationships at the polyline endpoints
    if ( P_first_halfedge == next_marked_halfedge_around_target_vertex(P_previous, border_halfedges) )
    {
      decorator.set_prev(O_first_halfedge, O_previous);
      O_previous->HBase::set_next(O_first_halfedge);
      decorator.set_prev(O_previous->opposite(), O_first_halfedge->opposite());
      O_first_halfedge->opposite()->HBase::set_next(O_previous->opposite());
    }
    else{
      CGAL_assertion("!OUCH see OUCHTAG");
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
    CGAL_assertion( Qhedge_to_Phedge.find(qhedge)!=Qhedge_to_Phedge.end() );
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
    const std::set<Vertex_handle>& interior_vertices_,
    typename Base::Hedge_map& Qhedge_to_Phedge_,
    PolyhedronPointPMap ppmap_
  ) : Base( Qhedge_to_Phedge_ )
    ,facets(facets_)
    ,interior_halfedges(interior_halfedges_)
    ,interior_vertices(interior_vertices_)
    ,ppmap(ppmap_)
  {}

  void operator()(HDS& hds)
  {
    HalfedgeDS_decorator<HDS> decorator(hds);
    std::vector<Halfedge_handle> hedge_to_update_vertices;
    typedef std::pair<Halfedge_handle, Halfedge_handle> Hedge_pair;
    BOOST_FOREACH(Hedge_pair p, Qhedge_to_Phedge)
    {
      hedge_to_update_vertices.push_back( p.second );
      hedge_to_update_vertices.push_back( p.second->opposite() ); // needed for non-closed polyline
    }

    //insert halfedges
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges){
      Halfedge_handle new_h = hds.edges_push_back( *h );
      Qhedge_to_Phedge[ h ] = new_h;

      //create a copy of interior vertices only once
      if (  h->vertex()->halfedge()==h &&
            interior_vertices.find(h->vertex()) != interior_vertices.end() )
      {
        Vertex_handle v =
          decorator.vertices_push_back( *(h->vertex()) );
        put(ppmap, v, get(ppmap, h->vertex() ) );
        decorator.set_vertex_halfedge(v, new_h);
        decorator.set_vertex(new_h, v);
        hedge_to_update_vertices.push_back( new_h );
      }
      if (  h->opposite()->vertex()->halfedge()==h->opposite() &&
            interior_vertices.find(h->opposite()->vertex()) != interior_vertices.end() )
      {
        Vertex_handle v =
          decorator.vertices_push_back( *(h->opposite()->vertex()) );
        put(ppmap, v, get(ppmap, h->opposite()->vertex() ) );
        decorator.set_vertex_halfedge(v, new_h->opposite());
        decorator.set_vertex(new_h->opposite(), v);
        hedge_to_update_vertices.push_back( new_h->opposite() );
      }
    }

    //create facets, connect halfedges, and create vertices
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

    //update incident vertices
    BOOST_FOREACH(Halfedge_handle h, hedge_to_update_vertices)
    {
      Vertex_handle v=h->vertex();
      Halfedge_handle next_around_vertex=h;
      do{
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
        #warning OUCHTAG Handle this properly: this appends only when having a polyhedron featuring only polylines \
                 during the import of patches from P, at end point of polylines that are not cycles
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
        CGAL_assertion (next_around_vertex->next() != Halfedge_handle());
        next_around_vertex=next_around_vertex->next()->opposite();
        decorator.set_vertex(next_around_vertex, v);
      }while(h != next_around_vertex);
    }
  }
};

template <class Polyhedron>
struct Intersection_polylines{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  const std::vector<Halfedge_handle>& P;
  const std::vector<Halfedge_handle>& Q;
  const std::vector<int>& lengths;
  Intersection_polylines(
    const std::vector<Halfedge_handle>& P_polylines,
    const std::vector<Halfedge_handle>& Q_polylines,
    const std::vector<int>& lengths_
  ) : P( P_polylines )
    , Q( Q_polylines )
    , lengths( lengths_ )
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
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning I am rather sure we can avoid this and set prev above see what I did below
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    //check if there are border hedges, if which case non-interior halfedges
    //sharing a source with a border hedge does not have a valid previous.
    BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
      if(h->is_border_edge())
      {
        if (h->is_border()) h=h->opposite();
        if (    interior_vertices.find(h->vertex())
             == interior_vertices.end() )
        {
          CGAL_assertion( h->next()->is_border() );//we marked it above!
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

    // put the halfedges on the boundary of the patch on the boundary of the polyhedron
    Face_handle border;
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      new_patch_border.push_back(hds.edges_push_back(*(h->opposite())));
      decorator.set_face(h, border);
      decorator.set_face(new_patch_border.back(), border);
    }

    // set next/prev relationship of border halfedges
    BOOST_FOREACH(Halfedge_handle h, patch_border_halfedges)
    {
      Halfedge_handle next=h->next();
      while(!next->is_border())
        next=next->opposite()->next();
      h->HBase::set_next(next);
      decorator.set_prev(next,h);

      // setting prev is only needed in case the polyhedron has a boundary
      // and the intersection polyline intersects its boundary
      if ( !h->prev()->is_border() || h->prev()->next()!=h ){
        Halfedge_handle prev=h->prev();
        do{
          prev=prev->opposite()->prev();
        } while( !prev->is_border() );
        prev->HBase::set_next(h);
        decorator.set_prev(h, prev);
      }

      CGAL_assertion( h->prev()->is_border() );

      decorator.set_vertex_halfedge(h->vertex(),h);
    }

    BOOST_FOREACH(Halfedge_handle h, new_patch_border)
    {
      Halfedge_handle next=h->next();
      while(!next->is_border())
        next=next->opposite()->next();
      h->HBase::set_next(next);
      decorator.set_prev(next,h);

      // setting prev is only needed in case the polyhedron has a boundary
      // and the intersection polyline intersects its boundary
      if ( !h->prev()->is_border() || h->prev()->next()!=h ){
        Halfedge_handle prev=h->prev();
        do {
          prev=prev->opposite()->prev();
        } while( !prev->is_border() );
        prev->HBase::set_next(h);
        decorator.set_prev(h, prev);
      }

      CGAL_assertion( h->prev()->is_border() );
    }

    /*
    if (clean_up_patch)
    {
      //now remove the simplices
      BOOST_FOREACH(Halfedge_handle h, interior_halfedges)
        hds.edges_erase(h);
      BOOST_FOREACH(Face_handle f, facets)
        hds.faces_erase(f);
      BOOST_FOREACH(Vertex_handle v, interior_vertices)
        hds.vertices_erase(v);
    }
    */
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
    if ( patch_ids[ facet_indices[fit] ]==patch_id )
    {
      facets.push_back( fit );
      Halfedge_handle hedges[]={fit->halfedge(), fit->halfedge()->next(), fit->halfedge()->prev()};
      for (int i=0;i<3;++i)
      {
        if ( is_marked_edge.find(hedges[i])==is_marked_edge.end() )
        {
          if ( hedges[i] < hedges[i]->opposite() )
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
};

} //end of namespace internal_IOP

namespace Corefinement
{

template <class Polyhedron,
          class Facet_id_pmap_=Default,
          class Kernel_=Default,
          class PolyhedronPointPMap_=Default >
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
    typename boost::property_map<Polyhedron, boost::edge_index_t>::type
  >::type                                                         Facet_id_pmap;
public:
//Boolean operation indices
  enum {P_UNION_Q = 0, P_INTER_Q, P_MINUS_Q, Q_MINUS_P };
private:
//Data members
  Polyhedron *P_ptr, *Q_ptr;
  cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output;
  Facet_id_pmap P_facet_id_pmap, Q_facet_id_pmap;
  PolyhedronPointPMap ppmap;
  // bitset containing information about operations that cannot be
  // performed because of non-manifoldness or that is ambiguous
  // 0 = P+Q
  // 1 = P inter Q
  // 2 = P - Q
  // 3 = Q  - P
  std::bitset<4> impossible_operation;


//Polyhedron typedefs
  typedef typename Polyhedron::Halfedge_const_handle      Halfedge_const_handle;
  typedef typename Polyhedron::Halfedge_handle                  Halfedge_handle;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
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
    try{
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
      // P_ptr->is_valid() might not be true at this point. For example if two
      // incident patches are removed, there is a cycle of edges with not facets
      // on both side (one of the two patches is a coplanar patch).
    }
  }

  void disconnect_patches_from_polyhedra(
    Polyhedron* P_ptr,
    const boost::dynamic_bitset<>& patches_to_remove,
    Patch_container& patches_of_P
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
      patch.patch_border_halfedges.swap(new_patch_border);
      // P_ptr->is_valid() might not be true at this point. For example if two
      // incident patches are removed, there is a cycle of edges with not facets
      // on both side (one of the two patches is a coplanar patch).
    }
  }

  template <bool reverse_patch_orientation >
  void append_Q_patches_to_P(
    Polyhedron* P_ptr,
    Polyhedron* Q_ptr,
    const boost::dynamic_bitset<>& patches_to_append,
    const std::vector<std::size_t>& Q_patch_ids,
    Facet_id_pmap Q_facet_id_pmap,
    const std::map< Halfedge_const_handle,
                      std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    std::map< Halfedge_handle,
            Halfedge_handle,
            internal_IOP::Compare_unik_address<Polyhedron>
          >& Qhedge_to_Phedge
  ){
    if ( patches_to_append.any() )
    {
      std::vector<Facet_handle> facets;
      std::set<Vertex_handle> interior_vertices;
      std::vector<Halfedge_handle> interior_halfedges;
      std::vector<Halfedge_handle> patch_border_halfedges;

      for (std::size_t i=patches_to_append.find_first();
                       i < patches_to_append.npos; i = patches_to_append.find_next(i))
      {
        internal_IOP::extract_patch_simplices(i, *Q_ptr, Q_facet_id_pmap, Q_patch_ids, facets, interior_vertices, interior_halfedges, patch_border_halfedges, border_halfedges);
      }

      //CGAL_assertion( Qhedge_to_Phedge.size()==patch_border_halfedges.size() ); // this is not true in case of coplanar patches
      internal_IOP::Surface_extension_by_patch_appending<Polyhedron,PolyhedronPointPMap, reverse_patch_orientation>
        modifier(facets, interior_halfedges, interior_vertices, Qhedge_to_Phedge, ppmap);
      P_ptr->delegate(modifier);
      //CGAL_assertion( P_ptr->is_valid() );
    }
  }

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
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges)
  {
    internal_IOP::Import_polyline<Polyhedron, PolyhedronPointPMap>
      modifier( vertex_map, P_to_O_hedge, Q_to_O_hedge, nb_segments,
                P_first_polyline_hedge, Q_first_polyline_hedge,
                border_halfedges, ppmap);

    O.delegate( modifier );
  }

  void compute_difference_inplace(
    Polyhedron* P_ptr,
    Polyhedron* Q_ptr,
    boost::dynamic_bitset<> is_patch_inside_Q,
    boost::dynamic_bitset<> is_patch_inside_P,
    const boost::dynamic_bitset<>& coplanar_patches_of_P,
    const boost::dynamic_bitset<>& coplanar_patches_of_Q,
    const boost::dynamic_bitset<>& coplanar_patches_of_P_for_union_and_intersection,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    const Intersection_polylines& polylines
  ){
      std::map< Halfedge_handle,
                Halfedge_handle,
                internal_IOP::Compare_unik_address<Polyhedron>
              > Qhedge_to_Phedge;
      //maps patch border halfedge from Q to halfedge from P
      std::size_t nb_polylines = polylines.lengths.size();
      for( std::size_t i=0; i<nb_polylines; ++i)
      {
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

      if ( coplanar_patches_of_P.any() )
        is_patch_inside_Q|=(coplanar_patches_of_P & coplanar_patches_of_P_for_union_and_intersection);

      //clean up patches inside Q
      remove_patches_from_polyhedra(P_ptr, is_patch_inside_Q, patches_of_P);

      //we import patches from Q
      boost::dynamic_bitset<> is_patch_outside_P=is_patch_inside_P;
      is_patch_outside_P.flip();
      is_patch_inside_P-=coplanar_patches_of_Q;
      append_Q_patches_to_P<true>(P_ptr, Q_ptr, is_patch_inside_P, patches_of_Q.patch_ids, patches_of_Q.facet_id_pmap, patches_of_Q.border_halfedges, Qhedge_to_Phedge);

      /// in case the result is empty, there will be no facets in the polyhedron
      /// but maybe marked halfedges
      if ( P_ptr->facets_begin()==P_ptr->facets_end() )
        P_ptr->clear();
  }
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
  #warning factorize with the code of compute_difference_inplace
#endif
  void compute_difference_inplace_delay_removal(
    Polyhedron* P_ptr,
    Polyhedron* Q_ptr,
    boost::dynamic_bitset<> is_patch_inside_Q,
    boost::dynamic_bitset<> is_patch_inside_P,
    const boost::dynamic_bitset<>& coplanar_patches_of_P,
    const boost::dynamic_bitset<>& coplanar_patches_of_Q,
    const boost::dynamic_bitset<>& coplanar_patches_of_P_for_union_and_intersection,
    Patch_container& patches_of_P,
    Patch_container& patches_of_Q,
    const Intersection_polylines& polylines
  ){
      std::map< Halfedge_handle,
                Halfedge_handle,
                internal_IOP::Compare_unik_address<Polyhedron>
              > Qhedge_to_Phedge;
      //maps patch border halfedge from Q to halfedge from P
      std::size_t nb_polylines = polylines.lengths.size();
      for( std::size_t i=0; i<nb_polylines; ++i)
      {
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

      if ( coplanar_patches_of_P.any() )
        is_patch_inside_Q|=(coplanar_patches_of_P & coplanar_patches_of_P_for_union_and_intersection);

      //clean up patches inside Q
      disconnect_patches_from_polyhedra(P_ptr, is_patch_inside_Q, patches_of_P);

      //we import patches from Q
      boost::dynamic_bitset<> is_patch_outside_P=is_patch_inside_P;
      is_patch_outside_P.flip();
      is_patch_inside_P-=coplanar_patches_of_Q;
      append_Q_patches_to_P<true>(P_ptr, Q_ptr, is_patch_inside_P, patches_of_Q.patch_ids, patches_of_Q.facet_id_pmap, patches_of_Q.border_halfedges, Qhedge_to_Phedge);

      /// in case the result is empty, there will be no facets in the polyhedron
      /// but maybe marked halfedges
      if ( P_ptr->facets_begin()==P_ptr->facets_end() )
        P_ptr->clear();
  }

  void fill_new_polyhedron(
    Polyhedron& O, // output
    const boost::dynamic_bitset<>& patches_of_P,
    const boost::dynamic_bitset<>& patches_of_Q,
    const std::vector<std::size_t>& P_patch_ids,
    const std::vector<std::size_t>& Q_patch_ids,
    Facet_id_pmap P_facet_id_pmap,
    Facet_id_pmap Q_facet_id_pmap,
    const Intersection_polylines& polylines,
    const std::map< Halfedge_const_handle,
                    std::pair<int,int>,Cmp_unik_ad >& border_halfedges
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
        import_polyline(O, polylines.P[i], polylines.Q[i],
                        polylines.lengths[i],
                        P_to_O_hedge, Q_to_O_hedge,
                        P_to_O_vertex, border_halfedges);

    //import patches of P
    append_Q_patches_to_P<false>(&O, P_ptr, patches_of_P, P_patch_ids, P_facet_id_pmap, border_halfedges, P_to_O_hedge);

    //import patches from Q
    append_Q_patches_to_P<false>(&O, Q_ptr, patches_of_Q, Q_patch_ids, Q_facet_id_pmap, border_halfedges, Q_to_O_hedge);

    /// in case the result is empty, there will be no facets in the polyhedron
    /// but maybe marked halfedges
    if ( O.facets_begin()==O.facets_end() )
      O.clear();
  }

public:
  static const bool do_compute_intersection_of_coplanar_facets = false;

  Polyhedra_output_builder(
    Polyhedron& P,
    Polyhedron& Q,
    cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output_,
    Facet_id_pmap P_facet_id_pmap_,
    Facet_id_pmap Q_facet_id_pmap_,
    PolyhedronPointPMap point_pmap = PolyhedronPointPMap()
  ) : P_ptr(&P)
    , Q_ptr(&Q)
    , desired_output( desired_output_ )
    , P_facet_id_pmap(P_facet_id_pmap_)
    , Q_facet_id_pmap(Q_facet_id_pmap_)
    , ppmap(point_pmap) {}

    Polyhedra_output_builder(
    Polyhedron& P,
    Polyhedron& Q,
    cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output_,
    PolyhedronPointPMap point_pmap = PolyhedronPointPMap()
  ) : P_ptr(&P)
    , Q_ptr(&Q)
    , desired_output( desired_output_ )
    , P_facet_id_pmap( boost::get(boost::edge_index, P) )
    , Q_facet_id_pmap( boost::get(boost::edge_index, Q) )
    , ppmap(point_pmap) {}

  bool union_valid()        const { return !impossible_operation[P_UNION_Q]; }
  bool intersection_valid() const { return !impossible_operation[P_INTER_Q]; }
  bool P_minus_Q_valid()    const { return !impossible_operation[P_MINUS_Q]; }
  bool Q_minus_P_valid()    const { return !impossible_operation[Q_MINUS_P]; }

  template <class Nodes_vector>
  void operator()(
    const std::map<Halfedge_const_handle,
                   std::pair<int,int>,Cmp_unik_ad >& border_halfedges,
    const Nodes_vector& nodes,
    const An_edge_per_polyline_map& an_edge_per_polyline,
    const Poly_to_map_node& /* polyhedron_to_map_node_to_polyhedron_vertex */)
  {

    cpp11::array< boost::optional<Polyhedron*>, 4 > results;
    // if (results[0]), it will hold P + Q
    // if (results[1]), it will hold P inter Q
    // if (results[2]), it will hold P - Q
    // if (results[3]), it will hold Q - P

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
    init_facet_indices(*P_ptr, P_facet_id_pmap);
    std::size_t P_nb_cc=
      mark_connected_components_v2(*P_ptr,
                                   is_not_marked,
                                   P_facet_id_pmap,
                                   P_patch_ids,
                                   P_patch_sizes);

    init_facet_indices(*Q_ptr, Q_facet_id_pmap);
    std::size_t Q_nb_cc=
      mark_connected_components_v2(*Q_ptr,
                                   is_not_marked,
                                   Q_facet_id_pmap,
                                   Q_patch_ids,
                                   Q_patch_sizes);

    CGAL_assertion ( P_nb_cc==P_patch_sizes.size() );
    CGAL_assertion ( Q_nb_cc==Q_patch_sizes.size() );

    // (2-a) Use the orientation around an edge to classify a patch
    boost::dynamic_bitset<> is_patch_inside_Q(P_nb_cc, false);
    boost::dynamic_bitset<> is_patch_inside_P(Q_nb_cc, false);
    boost::dynamic_bitset<> patch_status_not_set_P(P_nb_cc,true);
    boost::dynamic_bitset<> patch_status_not_set_Q(Q_nb_cc,true);
    boost::dynamic_bitset<> coplanar_patches_of_P(P_nb_cc,false);
    boost::dynamic_bitset<> coplanar_patches_of_Q(Q_nb_cc,false);
    boost::dynamic_bitset<> coplanar_patches_of_P_for_union_and_intersection(P_nb_cc,false);
    boost::dynamic_bitset<> coplanar_patches_of_Q_for_difference(Q_nb_cc,false);

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
            //Union allowed everything else is non_manifold
            impossible_operation.set();
            impossible_operation.reset(0);

            std::size_t patch_id_P =
              P_patch_ids[ P_facet_id_pmap[ first_hedge->is_border()
                                            ? first_hedge->opposite()->facet()
                                            : first_hedge->facet() ] ];
            std::size_t patch_id_Q =
              Q_patch_ids[ Q_facet_id_pmap[ second_hedge->is_border()
                                            ? second_hedge->opposite()->facet()
                                            : second_hedge->facet() ] ];
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
          int index_p1=node_index_of_incident_vertex(first_hedge->opposite()->next(),border_halfedges);
          int index_p2=node_index_of_incident_vertex(first_hedge->next(),border_halfedges);
          int index_q1=node_index_of_incident_vertex(second_hedge->opposite()->next(),border_halfedges);
          int index_q2=node_index_of_incident_vertex(second_hedge->next(),border_halfedges);

          std::size_t patch_id_p1=P_patch_ids[ P_facet_id_pmap[first_hedge->opposite()->facet()] ];
          std::size_t patch_id_p2=P_patch_ids[ P_facet_id_pmap[first_hedge->facet()] ];
          std::size_t patch_id_q1=Q_patch_ids[ Q_facet_id_pmap[second_hedge->opposite()->facet()] ];
          std::size_t patch_id_q2=Q_patch_ids[ Q_facet_id_pmap[second_hedge->facet()] ];

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

            bool Q2_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);
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
              coplanar_patches_of_Q_for_difference.set(patch_id_q2);
              bool Q1_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
              if ( Q1_is_between_P1P2 ) is_patch_inside_P.set(patch_id_q1); //case 3
              // else case 4
              continue;
            }
            else
            {
              if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p2,index_q1,P2,Q1,nodes) ) //P2==Q1
              {
                coplanar_patches_of_P.set(patch_id_p2);
                coplanar_patches_of_Q.set(patch_id_q1);
                coplanar_patches_of_Q_for_difference.set(patch_id_q1);
                bool Q2_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);
                if ( Q2_is_between_P1P2 ) is_patch_inside_P.set(patch_id_q2); //case 5
                // else case 6
                continue;
              }
              else{
                if ( are_triangles_coplanar_same_side_filtered(indices.first,indices.second,index_p2,index_q2,P2,Q2,nodes) ) //P2==Q2
                {
                  coplanar_patches_of_P.set(patch_id_p2);
                  coplanar_patches_of_Q.set(patch_id_q2);
                  coplanar_patches_of_P_for_union_and_intersection.set(patch_id_p2);
                  bool Q1_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
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

          CGAL_assertion(get(ppmap,P1) !=get(ppmap,Q1) && get(ppmap,P1)!=get(ppmap,Q2) && get(ppmap,P2) !=get(ppmap,Q1) && get(ppmap,P2)!=get(ppmap,Q2));

          bool Q1_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q1,P1,P2,Q1,nodes,ppmap);
          bool Q2_is_between_P1P2 = sorted_around_edge_filtered(indices.first,indices.second,index_p1,index_p2,index_q2,P1,P2,Q2,nodes,ppmap);

          if ( Q1_is_between_P1P2 ){
            is_patch_inside_P.set(patch_id_q1);
            if( Q2_is_between_P1P2 )
            {
              is_patch_inside_P.set(patch_id_q2);
              bool P1_is_between_Q1Q2 = sorted_around_edge_filtered(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes,ppmap);
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
              is_patch_inside_P.set(patch_id_q2);
              is_patch_inside_Q.set(patch_id_p1);
            }
            else
            {
              bool P1_is_between_Q1Q2 = sorted_around_edge_filtered(indices.first,indices.second,index_q1,index_q2,index_p1,Q1,Q2,P1,nodes,ppmap);
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
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    typedef Triangle_accessor_with_ppmap_3<Polyhedron,
                                           PolyhedronPointPMap> T3_accessor;
    T3_accessor accessor(ppmap);
    typedef Point_inside_polyhedron_3<Polyhedron, Kernel, T3_accessor> Inside_poly_test;
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning we must do the following only if we need them, that is we take the Union or P-Q or Q-P (only one is needed then)\
             for example this is useless if we consider the intersection
    #warning this should not be done if we have surfaces with boundaries!!! ask the user a flag?

    #warning Maybe we should forbid having several CC??
    #warning in case of two identical meshes, we will end up directly here, we need to a have a special handling
#endif
    if ( an_edge_per_polyline.begin()==an_edge_per_polyline.end() )
    {
      // the models are either identical (two indentical meshes or same surface, like two cubes meshed differently) or disjoint
      // \todo build the AABB-tree
      // call any intersected primitive with a point, if no intersection then disjoint, otherwise same model
    }

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning stop using next_marked_halfedge_around_target_vertex and create lists of halfedges instead?

    #warning this test is not robust if we have several CC and one has intersection and the other do not (thus forbidding several CC)
#endif

    if ( patch_status_not_set_P.any() )
    {
      Inside_poly_test inside_Q(*Q_ptr, accessor);
      for (typename Polyhedron::Face_iterator fit=P_ptr->facets_begin(),
                                              fit_end=P_ptr->facets_end();
                                              fit!=fit_end; ++fit)
      {
        std::size_t patch_id=P_patch_ids[ P_facet_id_pmap[fit] ];
        if ( patch_status_not_set_P.test( patch_id ) )
        {
          patch_status_not_set_P.reset( patch_id );
          if ( inside_Q( get(ppmap,fit->halfedge()->vertex()) ) == ON_BOUNDED_SIDE )
            is_patch_inside_Q.set(patch_id);

          if ( patch_status_not_set_P.none() ) break;
        }
      }
    }

    if ( patch_status_not_set_Q.any() )
    {
      Inside_poly_test inside_P(*P_ptr, accessor);
      for (typename Polyhedron::Face_iterator fit=Q_ptr->facets_begin(),
                                              fit_end=Q_ptr->facets_end();
                                              fit!=fit_end; ++fit)
      {
        std::size_t patch_id=Q_patch_ids[ Q_facet_id_pmap[fit] ];
        if ( patch_status_not_set_Q.test( patch_id ) )
        {
          patch_status_not_set_Q.reset( patch_id );
          if ( inside_P( get(ppmap,fit->halfedge()->vertex()) ) == ON_BOUNDED_SIDE )
            is_patch_inside_P.set(patch_id);

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
    std::cout << "coplanar_patches_of_Q_for_difference " << coplanar_patches_of_Q_for_difference << "\n";
    std::cout << "Size of patches of P: ";
    std::copy(P_patch_sizes.begin(), P_patch_sizes.end(), std::ostream_iterator<std::size_t>(std::cout," ") );
    std::cout << "\n";
    std::cout << "Size of patches of Q: ";
    std::copy(Q_patch_sizes.begin(), Q_patch_sizes.end(), std::ostream_iterator<std::size_t>(std::cout," ") );
    std::cout << "\n";
    #endif

    //For now I compute the union and the intersection is new polyhedra
    //and P and Q are updated to respectively contain P-Q and Q-P

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
      polyline_lengths.push_back(polyline_info.second-1);
    }
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning Factorize patch extraction
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG

    //store the patch description in a container to avoid recomputing it several times
    Patch_container patches_of_P( P_ptr, P_patch_ids, P_facet_id_pmap, border_halfedges, P_nb_cc),
                    patches_of_Q( Q_ptr, Q_patch_ids, Q_facet_id_pmap, border_halfedges, Q_nb_cc);

    /// compute the union
    if ( !impossible_operation.test(P_UNION_Q) && desired_output[P_UNION_Q] )
    {
      //tmp: compute the union in a new polyhedra
      Polyhedron* union_ptr = *desired_output[P_UNION_Q];
      CGAL_assertion(P_ptr!=union_ptr && Q_ptr!=union_ptr); //tmp

      //define patches to import from P
      boost::dynamic_bitset<> is_patch_outside_Q=is_patch_inside_Q;
      if ( coplanar_patches_of_P.any() )
        is_patch_outside_Q|=(coplanar_patches_of_P - coplanar_patches_of_P_for_union_and_intersection);
      is_patch_outside_Q.flip();

      //define patches to import from Q
      boost::dynamic_bitset<> is_patch_outside_P=is_patch_inside_P;
      is_patch_outside_P.flip();
      is_patch_outside_P-=coplanar_patches_of_Q;

      fill_new_polyhedron(
        *union_ptr,
        is_patch_outside_Q, is_patch_outside_P,
        P_patch_ids, Q_patch_ids,
        P_facet_id_pmap, Q_facet_id_pmap,
        Intersection_polylines( P_polylines,
                                Q_polylines,
                                polyline_lengths),
        border_halfedges
      );

      // std::ofstream output("P_union_Q.off");
      // output << Union;
    }

    /// compute the intersection
    if ( !impossible_operation.test(P_INTER_Q) && desired_output[P_INTER_Q] )
    {
      // tmp: compute the intersection in a new polyhedra
      Polyhedron* inter_ptr = *desired_output[P_INTER_Q];
      CGAL_assertion(P_ptr!=inter_ptr && Q_ptr!=inter_ptr); //tmp

      fill_new_polyhedron(
        *inter_ptr,
        coplanar_patches_of_P.any() ? is_patch_inside_Q
                                    : is_patch_inside_Q |
                                      (coplanar_patches_of_P - coplanar_patches_of_P_for_union_and_intersection),
        is_patch_inside_P,
        P_patch_ids, Q_patch_ids,
        P_facet_id_pmap, Q_facet_id_pmap,
        Intersection_polylines( P_polylines,
                                Q_polylines,
                                polyline_lengths),
        border_halfedges
      );

      // std::ofstream output("P_inter_Q.off");
      // output << Inter;
    }

    /// tmp: update P to contain P-Q
    if ( !impossible_operation.test(P_MINUS_Q) && desired_output[P_MINUS_Q] )
    {
      CGAL_assertion( *desired_output[P_MINUS_Q] == P_ptr ); //tmp

      compute_difference_inplace(
      //~ compute_difference_inplace_delay_removal(
        P_ptr, Q_ptr,
        is_patch_inside_Q, is_patch_inside_P,
        coplanar_patches_of_P,
        coplanar_patches_of_Q,
        coplanar_patches_of_P_for_union_and_intersection,
        patches_of_P, patches_of_Q,
        Intersection_polylines(P_polylines, Q_polylines, polyline_lengths)
      );
      // std::ofstream output("P_minus_Q.off");
      // output << *P_ptr;
    }

    /// tmp: update Q to contain Q-P
    if ( !impossible_operation.test(Q_MINUS_P) && desired_output[Q_MINUS_P] )
    {
#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
      #warning We must not clean patches from P, simply disconnect them and only remove them at the end.
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG

      CGAL_assertion( *desired_output[Q_MINUS_P] == Q_ptr ); //tmp

      compute_difference_inplace( Q_ptr, P_ptr,
                                  is_patch_inside_P, is_patch_inside_Q,
                                  coplanar_patches_of_Q,
                                  coplanar_patches_of_P,
                                  ~coplanar_patches_of_Q_for_difference & coplanar_patches_of_Q,
                                  patches_of_Q, patches_of_P,
                                  Intersection_polylines(Q_polylines, P_polylines, polyline_lengths)
      );
      // std::ofstream output("Q_minus_P.off");
      // output << *Q_ptr;
    }
  }
};

} } // end of namespace CGAL::Corefinement

#endif // CGAL_INTERNAL_COREFINEMENT_COMBINATORIAL_MAP_OUTPUT_BUILDER_H