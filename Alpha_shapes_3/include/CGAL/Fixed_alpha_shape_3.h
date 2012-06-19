// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the so
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_FIXED_ALPHA_SHAPE_3_H
#define CGAL_FIXED_ALPHA_SHAPE_3_H

#include <CGAL/basic.h>

#include <set>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <queue>
#include <boost/next_prior.hpp>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/iterator.h>
#ifdef CGAL_USE_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>  // TBC
#endif

#include <CGAL/internal/Classification_type.h>

#include <CGAL/Triangulation_3.h>


namespace CGAL{

namespace internal{
  //small utility to select the correct predicate in the weighted case
  template <class One_alpha,class Weighted_tag>
  struct Simplex_classif_predicate;

  template <class One_alpha>
  struct Simplex_classif_predicate<One_alpha,Tag_false>{
    static
    typename One_alpha::Triangulation::Geom_traits::Compare_squared_radius_3
    predicate(const typename One_alpha::Triangulation& T){
      return T.geom_traits().compare_squared_radius_3_object();
    }
  };

  template <class One_alpha>
  struct Simplex_classif_predicate<One_alpha,Tag_true>{
    static
    typename One_alpha::Triangulation::Geom_traits::Compare_weighted_squared_radius_3
    predicate(const typename One_alpha::Triangulation& T){
      return T.geom_traits().compare_weighted_squared_radius_3_object();
    }
  };
  
  //small utility to insert hidden vertices after a removal (in non-weighted case do nothing)
  template <class One_alpha,class Weighted_tag=::CGAL::Tag_true>
  struct Hidden_inserter
  {
    typedef typename One_alpha::Triangulation Dt;
    
    template <class Vertex_remover>
    static void insert(One_alpha& one_alpha, Vertex_remover& remover)
    {
      typename One_alpha::Cell_handle c;
      for (typename Vertex_remover::Hidden_points_iterator
        hi = remover.hidden_points_begin();
        hi != remover.hidden_points_end(); ++hi) 
      {
        typename One_alpha::Vertex_handle hv = one_alpha.insert (*hi, c);
        if (hv != typename One_alpha::Vertex_handle()) c = hv->cell();
      }
    }
  };
  
  template <class One_alpha>
  struct Hidden_inserter<One_alpha,::CGAL::Tag_false>{
    typedef typename One_alpha::Triangulation Dt;
    template <class Vertex_remover>
    static void insert(const One_alpha&,const Vertex_remover&){}
  };
  
  
} //internal


template < class Dt >
class Fixed_alpha_shape_3 : public Dt
{
  // DEFINITION The class Fixed_alpha_shape_3<Dt> represents the 
  // alpha-shape for a set of points (or a set of weighted points)
  // for a given value of alpha. The alphashape is defined  through
  // the Delaunay tetrahedralization of the points
  // (or the Regular tetrahedralization in case of weighted points)
  // and depends on the value of a parameter called alpha.
  // The alpha_shape is the domain of a subcomplex of this triangulation
  // called the Alpha_complex. The alpha_complex includes any simplex
  // having  a circumscribing sphere (an orthogonal sphere
  // in case of weighted points) empty of other points
  // (or suborthogonal to other sites in case of weighted points)
  // with squared radius equal or less than alpha
 
  // In each k-dimensional simplex of the triangulation
  // for (k=0,1,2) 
  // can be classified as EXTERIOR, SINGULAR, REGULAR
  // or INTERIOR with respect to the alpha shape.
  // A $k$ simplex is REGULAR if it is on the boundary
  // of the alpha_complex and belongs to a $k+1$ simplex in the complex
  // and it is SINGULAR simplex if it is  a boundary simplex tht is not
  // included in a $k+1$ simplex of the complex.
  
  // Roughly, the Fixed_alpha_shape data structure computes and stores, 
  // for each simplex it classification type.
  

  //------------------------- TYPES ------------------------------------

public:
  typedef Dt                                        Triangulation;
  typedef typename Dt::Geom_traits                  Gt;
  typedef typename Dt::Triangulation_data_structure Tds;

  //Classification type: no longer an enum inside the class as Vertex and cell must know it
  typedef internal::Classification_type             Classification_type;
  static const Classification_type                  EXTERIOR = internal::EXTERIOR;
  static const Classification_type                  REGULAR  = internal::REGULAR;
  static const Classification_type                  INTERIOR = internal::INTERIOR;
  static const Classification_type                  SINGULAR = internal::SINGULAR;

  typedef typename Gt::FT                           Coord_type;
  typedef Coord_type                                NT;
  typedef Coord_type                                FT;

  typedef typename Gt::Point_3                      Point;
  
  typedef typename Dt::Cell_handle                  Cell_handle;
  typedef typename Dt::Vertex_handle                Vertex_handle;
  typedef typename Dt::Facet                        Facet;
  typedef typename Dt::Edge                         Edge;

  typedef typename Dt::Cell_circulator              Cell_circulator;
  typedef typename Dt::Facet_circulator             Facet_circulator;

  typedef typename Dt::Cell_iterator                Cell_iterator;
  typedef typename Dt::Facet_iterator               Facet_iterator;
  typedef typename Dt::Edge_iterator                Edge_iterator;
  typedef typename Dt::Vertex_iterator              Vertex_iterator;

  typedef typename Dt::Finite_cells_iterator        Finite_cells_iterator;
  typedef typename Dt::Finite_facets_iterator       Finite_facets_iterator;
  typedef typename Dt::Finite_edges_iterator        Finite_edges_iterator;
  typedef typename Dt::Finite_vertices_iterator     Finite_vertices_iterator;

  typedef typename Dt::Locate_type                  Locate_type;
  typedef typename Dt::Weighted_tag                 Weighted_tag;

  using Dt::dimension;
  using Dt::finite_facets_begin;
  using Dt::finite_facets_end;
  using Dt::finite_edges_begin;
  using Dt::finite_edges_end;
  using Dt::all_edges_begin;
  using Dt::all_edges_end;
  using Dt::finite_vertices_begin;
  using Dt::finite_vertices_end;
  using Dt::finite_cells_begin;
  using Dt::finite_cells_end;
  using Dt::VERTEX;
  using Dt::EDGE;
  using Dt::FACET;
  using Dt::CELL;
  using Dt::OUTSIDE_CONVEX_HULL;
  using Dt::OUTSIDE_AFFINE_HULL;
  using Dt::vertex_triple_index;
  using Dt::finite_adjacent_vertices;
  using Dt::find_conflicts;
  using Dt::is_edge;
  using Dt::incident_vertices;
  using Dt::incident_facets;
  using Dt::is_infinite;
  using Dt::is_Gabriel;
  using Dt::tds;

  typedef std::pair<Vertex_handle, Vertex_handle>          Vertex_handle_pair;
  typedef std::map<Vertex_handle_pair,Classification_type> Edge_status_map;

  //test if a simplex is exterior to the alpha-shape
  class Exterior_simplex_test{
    const Fixed_alpha_shape_3 * _as;
  public:
    Exterior_simplex_test() {}
    Exterior_simplex_test(const Fixed_alpha_shape_3 * as) {_as = as;}
    bool operator() ( const Finite_cells_iterator& fci) const {
      return _as->classify(fci) == EXTERIOR ;
    }
    bool operator() ( const Finite_vertices_iterator& fvi) const {
      return _as->classify(fvi) == EXTERIOR ;
    }
    bool operator() ( const Finite_facets_iterator& ffi) const {
      return _as->classify(*ffi) == EXTERIOR ;
    }
    bool operator() ( const Finite_edges_iterator& fei) const {
      return _as->classify(*fei) == EXTERIOR ;
    }
  };

  typedef Filter_iterator< Finite_vertices_iterator, Exterior_simplex_test> Alpha_shape_vertices_iterator;
  typedef Filter_iterator< Finite_facets_iterator, Exterior_simplex_test>   Alpha_shape_facets_iterator;
  typedef Filter_iterator< Finite_edges_iterator, Exterior_simplex_test>    Alpha_shape_edges_iterator;
  typedef Filter_iterator< Finite_cells_iterator, Exterior_simplex_test>    Alpha_shape_cells_iterator;
  
 
  Vertex_handle 
  insert(const Point& p,Cell_handle start=Cell_handle())
  {
    if (this->dimension() < 3){
      Vertex_handle new_v=Triangulation::insert(p,start);
      if (this->dimension() == 3) initialize_alpha();
      return new_v;
    }
    
    //Handle only case of dimension 3 of insert_in_conflict from Triangulation_3 class.
    typename Triangulation::Locate_type lt;
    int li, lj;
    Cell_handle c = this->locate(p, lt, li, lj, start);
    typename Triangulation::Conflict_tester_3 tester(p, this);
    if ((lt == VERTEX) &&
        (tester.compare_weight(c->vertex(li)->point(), p)==0) ) {
      return c->vertex(li);
    }
    // If the new point is not in conflict with its cell, it is hidden.
    if (!tester.test_initial_cell(c)) {
      this->hidden_point_visitor.hide_point(c,p);
      return Vertex_handle();
    }

    // Ok, we really insert the point now.
    // First, find the conflict region.
    std::vector<Cell_handle> cells;
    std::vector<Facet> facets_on_the_boundary_of_the_hole;
    
    cells.reserve(32);
    this->find_conflicts
      (c, tester, make_triple(std::back_inserter(facets_on_the_boundary_of_the_hole),
                              std::back_inserter(cells),
                              Emptyset_iterator()));

    Facet facet=*boost::prior(facets_on_the_boundary_of_the_hole.end());
    
    // Remember the points that are hidden by the conflicting cells,
    // as they will be deleted during the insertion.
    this->hidden_point_visitor.process_cells_in_conflict(cells.begin(), cells.end());


  //Before insertion:
    //recover edges on the boundary of the hole.
    std::set<Edge,Compare_edge> hole_boundary_edges;
    const int indices[3]={1,2,3};
    for (typename std::vector<Facet>::iterator 
      it=facets_on_the_boundary_of_the_hole.begin();
      it!=facets_on_the_boundary_of_the_hole.end();
      ++it)
    {
      Facet f=it->first->tds_data().is_in_conflict()?this->mirror_facet(*it):*it;
      CGAL_precondition(!f.first->tds_data().is_in_conflict());
      for (int i=0;i<3;++i){
        Edge edge(f.first,(indices[i]+f.second)%4,(indices[(i+1)%3]+f.second)%4);
        if (!this->is_infinite(edge))
          hole_boundary_edges.insert(edge);
      }
    }
    // Erase from edge_status_map, edges that will disappear:
    // they are not on the boudary of the hole
    std::set<Edge,Compare_edge> hole_edges;
    std::pair<typename std::set<Edge,Compare_edge>::iterator,bool> it_hedge_and_not_already_seen;
    for (typename std::vector<Cell_handle>::iterator it=cells.begin();it!=cells.end();++it){
      for (int i=0;i<3;++i)
        for (int k=i+1;k<4;++k)
        {
          Edge edge(*it,i,k);
          if (this->is_infinite(edge) || hole_boundary_edges.find(edge)!=hole_boundary_edges.end() ) continue;
          it_hedge_and_not_already_seen=hole_edges.insert(edge);
          if (!it_hedge_and_not_already_seen.second){
            edge_status_map.erase(make_vertex_handle_pair(*(it_hedge_and_not_already_seen.first))); //for infinite edges it does nothing
          }
        }
    }

  //Insertion using base triangulation
    Vertex_handle v = this->_insert_in_hole(p, cells.begin(), cells.end(),
                                            facet.first, facet.second);

  //After insertion
    //--set classification of new cells and facets that were on the boundary of the hole
    cells.clear();
    this->finite_incident_cells(v,std::back_inserter(cells));
    for (typename std::vector<Cell_handle>::iterator it=cells.begin();it!=cells.end();++it){
      set_cell_status(*it); //set cell status
      set_facet_classification_type( Facet( (*it),(*it)->index(v) ) ); //set facet status (incident to one cell of the hole)
    }
    //--set classification of new facets
    std::vector<Facet> facets;
    this->finite_incident_facets(v,std::back_inserter(facets));
    for (typename std::vector<Facet>::iterator fit=facets.begin();fit!=facets.end();++fit)
      set_facet_classification_type(*fit);
    //--init classif of new vertex
    v->set_classification_type(SINGULAR);
    std::vector<Edge> edges;
    this->finite_incident_edges(v,std::back_inserter(edges));
    //--set status of new edges + update status of new vertex
    for (typename std::vector<Edge>::iterator eit=edges.begin();eit!=edges.end();++eit){
      Classification_type status=compute_edge_status(eit->first, eit->second, eit->third);
      Vertex_handle_pair vhp = make_vertex_handle_pair( *eit );
      CGAL_precondition( edge_status_map.find(vhp)==edge_status_map.end() );
      edge_status_map.insert(std::make_pair(vhp, status));
      update_vertex_status(v,status);
    }
    //--set final status of new vertex
    Cell_handle tmp;
    int itmp1,itmp2;
    v->is_on_chull( this->is_edge(this->infinite_vertex(),v,tmp,itmp1,itmp2) );
    finalize_status_of_vertex(v);
    
    //--set classification of old edges
    for (typename std::set<Edge,Compare_edge>::iterator  
         eit=hole_boundary_edges.begin();eit!=hole_boundary_edges.end();++eit)
    {
      CGAL_precondition (!this->is_infinite(*eit));
      Classification_type status=compute_edge_status(eit->first, eit->second, eit->third);
      Vertex_handle_pair vhp = make_vertex_handle_pair( *eit );
      typename Edge_status_map::iterator it_status=edge_status_map.find(vhp);
      CGAL_precondition(it_status!=edge_status_map.end());
      it_status->second = status;
    }
    
    //--set status of old vertices + update is_on_chull
    //TODO: find a better way to do it : make an update
    std::vector<Vertex_handle> vertices;
    this->finite_adjacent_vertices(v,std::back_inserter(vertices));
    for (typename std::vector<Vertex_handle>::iterator vit=vertices.begin();vit!=vertices.end();++vit){
      if ( (*vit)->is_on_chull() )
      {
        (*vit)->is_on_chull( this->is_edge(this->infinite_vertex(),(*vit),tmp,itmp1,itmp2) );
      }
      set_vertex_status(*vit);
    }
    
    // Store the hidden points in their new cells.
    this->hidden_point_visitor.reinsert_vertices(v);
    return v;
  }
  
  void remove (Vertex_handle vertex_to_remove)
  {
    CGAL_precondition(vertex_to_remove!=Vertex_handle());
    CGAL_precondition(vertex_to_remove!=this->infinite_vertex());
    
    std::vector<Facet> link;
    std::vector<Vertex_handle> vertices_to_update;    
    std::map<Vertex_handle,Classification_type> old_classification;
    
    if (this->dimension() == 3)
    {
      //recover facet of the link: they are bounding
      //the hole made when removing the vertex
      std::vector<Cell_handle> incident_cells;
      this->finite_incident_cells(vertex_to_remove,std::back_inserter(incident_cells));
      for (typename std::vector<Cell_handle>::iterator it=
          incident_cells.begin();it!=incident_cells.end();++it)
      {
        int index=(*it)->index(vertex_to_remove);
        link.push_back( this->mirror_facet(Facet(*it,index)) );
        CGAL_assertion (!this->is_infinite(link.back()));
      }

      //get vertices that will need to be updated
      this->finite_adjacent_vertices(vertex_to_remove,std::back_inserter(vertices_to_update));
      
      //1-erase removed edges from edge_map
      //2-store old classification of vertices and set it to SINGULAR
      for(typename std::vector<Vertex_handle>::iterator it=vertices_to_update.begin();it!=vertices_to_update.end();++it){
        CGAL_precondition(edge_status_map.find(make_vertex_handle_pair(*it,vertex_to_remove)) != edge_status_map.end());
        edge_status_map.erase(make_vertex_handle_pair(*it,vertex_to_remove));
        old_classification.insert(std::make_pair(*it,(*it)->get_classification_type()));
        (*it)->set_classification_type(SINGULAR);
      }
    }
   
    //Do remove the vertex from the underlying triangulation
    Dt tmp;
    typename Dt:: template Vertex_remover<Dt> remover(tmp);
    typedef CGAL::Triangulation_3<typename Dt::Geom_traits,typename Dt::Triangulation_data_structure> Tr_Base;
    Tr_Base::remove(vertex_to_remove,remover);
    
    if (this->dimension()<3){
      edge_status_map.clear();
      return;
    }
    
    //recover new cells
    std::set<Cell_handle> new_cells;//cells in the hole
    std::set<Cell_handle> outer;//cells that are not in the hole
    std::queue<Cell_handle> to_visit;
    for (typename std::vector<Facet>::iterator it=link.begin();it!=link.end();++it){
      Facet mirror_facet=this->mirror_facet(*it);
      if ( !this->is_infinite(mirror_facet.first) )
      {
        to_visit.push( mirror_facet.first );
        outer.insert(it->first);
      }
    }
    
    while (!to_visit.empty()){
      Cell_handle cell=to_visit.front();
      to_visit.pop();
      if ( !new_cells.insert(cell).second ) continue; //already explored
      //explore neighbor cells
      for (int i=0;i<4;++i){
        Cell_handle candidate=cell->neighbor(i);
        if ( !this->is_infinite(candidate) && outer.find(candidate)==outer.end())
          to_visit.push(candidate);
      }
    }
    
    //set status of new cells
    for (typename std::set<Cell_handle>::iterator it=new_cells.begin();it!=new_cells.end();++it){
      CGAL_precondition(!this->is_infinite(*it));
      set_cell_status(*it);
    }
    
    //recover new facets + link facets
    //set vertex that are on the convex hull (those on a facet incident to infinite cell)
    while (!new_cells.empty())
    {
      Cell_handle cell=*new_cells.begin();
      new_cells.erase(new_cells.begin());
      for (int i=0;i<4;++i){
        if ( new_cells.find(cell->neighbor(i)) != new_cells.end() )
          link.push_back(Facet(cell,i));
        else{
          if ( this->is_infinite( cell->neighbor(i) ) ){
            link.push_back(Facet(cell,i));
            cell->vertex((i+1)%4)->is_on_chull(true);
            cell->vertex((i+2)%4)->is_on_chull(true);
            cell->vertex((i+3)%4)->is_on_chull(true);
          }
        }
      }
    }
    
    std::set<Edge,Compare_edge> new_edges;
    //1- set status of these facets
    //2- recover new edges + edges incident to link facets
    const int index[3]={1,2,3};
    for (typename std::vector<Facet>::iterator it=link.begin();it!=link.end();++it){
      set_facet_classification_type(*it);
      for (int i=0;i<3;++i){
        new_edges.insert(Edge(it->first,(it->second+index[i])%4,(it->second+index[(i+1)%3])%4));
      }
    }
   
    //set status of these edges
    for(typename std::set<Edge,Compare_edge>::iterator it=new_edges.begin();it!=new_edges.end();++it){
      Classification_type status=compute_edge_status(it->first, it->second, it->third);
      //cross links
      Vertex_handle_pair vhp = make_vertex_handle_pair( *it );
      edge_status_map[vhp]=status;
    
      //update status of incident vertices
      update_vertex_status(it->first->vertex(it->second),status);
      update_vertex_status(it->first->vertex(it->third),status);
    }
    
    //set final status of vertices
    for(typename std::vector<Vertex_handle>::iterator it=vertices_to_update.begin();it!=vertices_to_update.end();++it){
      //this one is working but should be more expensive
      //set_vertex_status(*it);      continue;
      
      Classification_type old_status=old_classification[*it];
      Classification_type status=(*it)->get_classification_type();
      
      CGAL_precondition( status!=SINGULAR ); // at least on edge is incident to that vertex
      
      if (status==INTERIOR){
        if (old_status!=INTERIOR || (*it)->is_on_chull())
          (*it)->set_classification_type(REGULAR);
        else
          (*it)->set_classification_type(INTERIOR);
        continue;
      }
      
      if (status==REGULAR) continue;

      if ( status==EXTERIOR && ( old_status==REGULAR || old_status==INTERIOR) ){ //if vertex was EXTERIOR or SINGULAR other edges are not in the alpha complex
        //check if former REGULAR vertex becomes EXTERIOR or SINGULAR
        std::list<Vertex_handle> incidentv;
        finite_adjacent_vertices(*it,back_inserter(incidentv));
        
        typename std::list<Vertex_handle>::iterator vvit=incidentv.begin();
        for( ; vvit != incidentv.end(); ++vvit) {
            //TODO: We take all edges -> WE SHOULD ONLY TAKE THOSE NOT IN THE HOLE
            Vertex_handle_pair vhp = make_vertex_handle_pair( *vvit, *it);
            CGAL_assertion(edge_status_map.find(vhp)!=edge_status_map.end());
            Classification_type status=edge_status_map[vhp];
            if (status!=EXTERIOR){
              (*it)->set_classification_type(REGULAR);
              break;
            }
        }
        
        if ( vvit != incidentv.end() ) continue;
      }
      
      if ( is_vertex_Gabriel((*it),Weighted_tag()) && is_gabriel_simplex_in_alpha_complex((*it)) )
        (*it)->set_classification_type(SINGULAR);
      else
        (*it)->set_classification_type(EXTERIOR);
    }
    
    //Insert possible hidden points
    internal::Hidden_inserter<Fixed_alpha_shape_3<Dt>,Weighted_tag>::insert(*this,remover);
  }
  
  
  
private:

  typedef internal::Simplex_classif_predicate<Fixed_alpha_shape_3<Dt>,Weighted_tag> Simplex_classif;


//data members
  NT _alpha; //the value of alpha
  Edge_status_map edge_status_map; //the map containing the status of edges


  //------------------------- CONSTRUCTORS ------------------------------
public:
  Fixed_alpha_shape_3(NT alpha=0):_alpha(alpha){} 

  Fixed_alpha_shape_3(Dt& dt, NT alpha = 0):_alpha(alpha){
    Dt::swap(dt);
    if (dimension() == 3) initialize_alpha();
  }
 
  // Introduces an alpha-shape `A' for the alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last
  template < class InputIterator >  
  Fixed_alpha_shape_3(const InputIterator& first,  
                const InputIterator& last,  
                const NT& alpha = 0): _alpha(alpha)
  {
    Dt::insert(first, last);
    if (dimension() == 3)	  initialize_alpha();
  }
 
 
  
private :

  template < class InputIterator >  
  int make_alpha_shape(const InputIterator& first, 
		       const InputIterator& last)
  {
    clear();
    int n = Dt::insert(first, last);
    if (dimension() == 3)	  initialize_alpha();
    return n;
  }

  //--------------------- INITIALIZATION OF PRIVATE MEMBERS -----------
  void initialize_status_of_cells();
  void initialize_status_of_facets();
  void initialize_status_of_edges();
  void initialize_status_of_vertices();
  void finalize_status_of_vertex(Vertex_handle);
  void initialize_alpha() {
    //identify vertices on the convex hull
    std::vector<Vertex_handle> chull;
    incident_vertices(this->infinite_vertex(),std::back_inserter(chull));
    for ( typename std::vector<Vertex_handle>::iterator it=chull.begin();it!=chull.end();++it)
      (*it)->is_on_chull(true);
    
    initialize_status_of_cells();
    initialize_status_of_facets();
    initialize_status_of_edges();
    initialize_status_of_vertices();
  }

private :
  // prevent default copy constructor and default assigment
  Fixed_alpha_shape_3(const Fixed_alpha_shape_3&);
  void operator=(const Fixed_alpha_shape_3&);  

  static
  Vertex_handle_pair
  make_vertex_handle_pair( Vertex_handle v1, Vertex_handle v2) {
    return v1 < v2 ? std::make_pair(v1,v2)
                   : std::make_pair(v2,v1);
  }

  static
  Vertex_handle_pair
  make_vertex_handle_pair( const Edge& e) {
    return make_vertex_handle_pair(e.first->vertex(e.second),e.first->vertex(e.third));
  }
  
  struct Compare_edge{
    bool operator()(const Edge& e1,const Edge& e2) const
    {
      return make_vertex_handle_pair(e1)<make_vertex_handle_pair(e2);
    }
  };
  
  bool is_vertex_Gabriel(Vertex_handle,Tag_false){return true;}
  template <class Vertexhandle>
  bool is_vertex_Gabriel(Vertexhandle v,Tag_true){return is_Gabriel(v);}
  void update_vertex_status(Vertex_handle v,Classification_type edge_status);
  void set_facet_classification_type(const Facet& f);
  void set_vertex_status(Vertex_handle v);
  inline void set_cell_status (Cell_handle c);
  Classification_type compute_edge_status( const Cell_handle&  c,int i,int j) const;
  
  //---------------------------------------------------------------------

public:
  // Returns the current alpha-value.
  const NT&  get_alpha() const
  {
    return _alpha;
  }
  
  void clear()
  {
    // clears the structure
    Dt::clear();
    edge_status_map.clear();
  }
  
  //--------------------- PREDICATES -----------------------------------

public:
  Classification_type  classify(const Cell_handle& s) const  {
    if (is_infinite(s)) return EXTERIOR;
    return s->get_classification_type(); 
  }
  Classification_type  classify(const Facet& f) const {
    if (is_infinite(f)) return EXTERIOR;
    return f.first->get_facet_classification_type(f.second); 
  }
  Classification_type  classify(const Edge& e) const  {
    if (is_infinite(e)) return EXTERIOR;
    return edge_status_map.find(make_vertex_handle_pair(e))->second; 
  }
  Classification_type  classify(const Vertex_handle& v) const {
    if (is_infinite(v)) return EXTERIOR;
    return v->get_classification_type(); 
  }


//------------------- GEOMETRIC PRIMITIVES ----------------------------
private:
  bool
  is_gabriel_simplex_in_alpha_complex (const Cell_handle& s) const{
    return
      Simplex_classif::predicate(*this)(
          s->vertex(0)->point(),
          s->vertex(1)->point(),
          s->vertex(2)->point(),
          s->vertex(3)->point(),
          get_alpha()
      ) !=POSITIVE;
  }

  bool
  is_gabriel_simplex_in_alpha_complex (const Cell_handle& s, const int& i) const{
    return
      Simplex_classif::predicate(*this)(
          s->vertex(vertex_triple_index(i,0))->point(),
          s->vertex(vertex_triple_index(i,1))->point(),
          s->vertex(vertex_triple_index(i,2))->point(),
          get_alpha()
      ) != POSITIVE;
  }

  bool
  is_gabriel_simplex_in_alpha_complex (const Facet& f) {
    return is_gabriel_simplex_in_alpha_complex(f.first, f.second);
  }

  bool
  is_gabriel_simplex_in_alpha_complex (const Cell_handle& s, const int& i, const int& j) const  {
    return 
      Simplex_classif::predicate(*this)(
          s->vertex(i)->point(),
          s->vertex(j)->point(),
          get_alpha()
      ) !=POSITIVE;
  }

  bool
  is_gabriel_simplex_in_alpha_complex (const Edge& e) const {
    return is_gabriel_simplex_in_alpha_complex(e.first,e.second,e.third);
  }

  bool
  is_gabriel_simplex_in_alpha_complex (const Vertex_handle& v) const {
    return
      Simplex_classif::predicate(*this)(
          v->point(),
          get_alpha()
      ) !=POSITIVE;
  }

  //---------------------------------------------------------------------
public:  
#ifdef CGAL_USE_GEOMVIEW
  void show_alpha_shape_faces(Geomview_stream &gv) const;
#endif


  //Iterators

  //---------------------------------------------------------------------
  Alpha_shape_vertices_iterator alpha_shape_vertices_begin() const {
    return CGAL::filter_iterator(finite_vertices_end(),Exterior_simplex_test(this),finite_vertices_begin());}
  Alpha_shape_vertices_iterator alpha_shape_vertices_end() const   {
    return CGAL::filter_iterator(finite_vertices_end(),Exterior_simplex_test(this));}
  //---------------------------------------------------------------------
  Alpha_shape_facets_iterator alpha_shape_facets_begin() const{
    return CGAL::filter_iterator(finite_facets_end(),Exterior_simplex_test(this),finite_facets_begin());}
  Alpha_shape_facets_iterator alpha_shape_facets_end() const {
    return CGAL::filter_iterator(finite_facets_end(),Exterior_simplex_test(this));}
  //---------------------------------------------------------------------
  Alpha_shape_edges_iterator alpha_shape_edges_begin() const{
    return CGAL::filter_iterator(finite_edges_end(),Exterior_simplex_test(this),finite_edges_begin());}
  Alpha_shape_edges_iterator alpha_shape_edges_end() const {
    return CGAL::filter_iterator(finite_edges_end(),Exterior_simplex_test(this));}
  //---------------------------------------------------------------------
  Alpha_shape_cells_iterator alpha_shape_cells_begin() const {
    return CGAL::filter_iterator(finite_cells_end(),Exterior_simplex_test(this),finite_cells_begin());
  }
  Alpha_shape_cells_iterator alpha_shape_cells_end() const {
    return CGAL::filter_iterator(finite_cells_end(),Exterior_simplex_test(this));
  }
  //---------------------------------------------------------------------
  // To extract simplices given a classification type
  template<class OutputIterator>
  OutputIterator get_alpha_shape_cells(OutputIterator it, 
				       Classification_type type) const
  {
    Finite_cells_iterator cit = finite_cells_begin();
    for( ; cit != finite_cells_end() ; ++cit){
      if (classify(cit) == type) *it++ = Cell_handle(cit);
    }
    return it;
  }

  template<class OutputIterator>
  OutputIterator get_alpha_shape_facets(OutputIterator it, 
					Classification_type type) const
  {
    Finite_facets_iterator fit = finite_facets_begin();
    for( ; fit != finite_facets_end() ; ++fit){
      if (classify(*fit) == type) *it++ = *fit;
    }
    return it;
  }

  template<class OutputIterator>
  OutputIterator get_alpha_shape_edges(OutputIterator it, 
				       Classification_type type) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for( ; eit != finite_edges_end() ; ++eit){
      if (classify(*eit) == type) *it++ = *eit;
    }
    return it;
  }

  template<class OutputIterator>
   OutputIterator get_alpha_shape_vertices(OutputIterator it, 
					   Classification_type type) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for( ; vit != finite_vertices_end() ; ++vit){
      if (classify(vit) == type) *it++ = Vertex_handle(vit);
    }
    return it;
  }
  
};


template < class Dt >
const typename Fixed_alpha_shape_3<Dt>::Classification_type Fixed_alpha_shape_3<Dt>::EXTERIOR;
template < class Dt >
const typename Fixed_alpha_shape_3<Dt>::Classification_type Fixed_alpha_shape_3<Dt>::REGULAR;
template < class Dt >
const typename Fixed_alpha_shape_3<Dt>::Classification_type Fixed_alpha_shape_3<Dt>::INTERIOR;
template < class Dt >
const typename Fixed_alpha_shape_3<Dt>::Classification_type Fixed_alpha_shape_3<Dt>::SINGULAR;

//---------------------------------------------------------------------
//--------------------- MEMBER FUNCTIONS-------------------------------
//---------------------------------------------------------------------


//--------------------- INITIALIZATION OF PRIVATE MEMBERS -------------

template <class Dt>
void Fixed_alpha_shape_3<Dt>::set_cell_status(Cell_handle c){
  Classification_type status=is_infinite(c) ? EXTERIOR:( is_gabriel_simplex_in_alpha_complex(c) ? INTERIOR : EXTERIOR );
  c->set_classification_type(status);  
}

template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::initialize_status_of_cells()
{ 
  Finite_cells_iterator cell_it, done = finite_cells_end();
  for( cell_it = finite_cells_begin(); cell_it != done; ++cell_it) {
    set_cell_status(cell_it);
  }
}


//---------------------------------------------------------------------

template < class Dt >
void
Fixed_alpha_shape_3<Dt>::
set_facet_classification_type( const Facet& f) {
  Cell_handle pCell = f.first;
  int i = f.second;
  Cell_handle pNeighbor = pCell->neighbor(i);
  int iNeigh = pNeighbor->index(pCell);

  unsigned nb_interior_cells=0;
  
  if(!is_infinite(pCell)){
    if (pCell->get_classification_type()==INTERIOR)
      ++nb_interior_cells;
  }
  
  if(!is_infinite(pNeighbor)){
    if (pNeighbor->get_classification_type()==INTERIOR)
      ++nb_interior_cells;
  }
  
  Classification_type status=EXTERIOR;
  switch (nb_interior_cells){
    case 2:
      status=INTERIOR;
    break;
    case 1:
      status=REGULAR;
    break;
    default:
    {
      if ( is_Gabriel(f) ){
        if ( is_gabriel_simplex_in_alpha_complex(f) )  status=SINGULAR;          
      }
    }
  }  
  pCell->set_facet_classification_type(i,status);
  pNeighbor->set_facet_classification_type(iNeigh,status);
}


template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::initialize_status_of_facets()
{
  for(Finite_facets_iterator fit = finite_facets_begin();
      fit != finite_facets_end(); ++fit)   
        set_facet_classification_type(*fit);
}



template < class Dt >
typename Fixed_alpha_shape_3<Dt>::Classification_type
Fixed_alpha_shape_3<Dt>::
compute_edge_status( const Cell_handle& c, int i, int j) const
{
  Facet_circulator fcirc, done;
  fcirc = incident_facets(c,i,j);
  done = fcirc;
  
  bool is_regular=false;
  bool is_interior=true;
  do{
    if (!is_infinite(*fcirc)){
      Classification_type status=(*fcirc).first->get_facet_classification_type((*fcirc).second);
      if (status!=INTERIOR) is_interior=false;
      if (status!=EXTERIOR)
        is_regular=true;
      else
        if (is_regular)
          break;
    }
  }while(++fcirc!=done);
  if (is_interior) return INTERIOR;
  if (is_regular)  return REGULAR;
  
  if ( is_Gabriel(c,i,j) ){
    if ( is_gabriel_simplex_in_alpha_complex(c,i,j) )  return SINGULAR;
  }  
  return EXTERIOR;
}

template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::initialize_status_of_edges()
{
  for (Finite_edges_iterator eit = finite_edges_begin(); 
       eit != finite_edges_end(); ++eit) 
  {
    Classification_type status=compute_edge_status(eit->first, eit->second, eit->third);
    //cross links
    Vertex_handle_pair vhp = make_vertex_handle_pair( *eit );
    edge_status_map.insert(std::make_pair(vhp, status));
  }
}


//this function is only to use for update (removal/update)
template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::set_vertex_status(Vertex_handle v){
  std::list<Vertex_handle> incidentv;
  finite_adjacent_vertices(v,back_inserter(incidentv));
  
  bool is_interior=true;
  bool is_regular=false;
  typename std::list<Vertex_handle>::iterator vvit=incidentv.begin();
  for( ; vvit != incidentv.end(); ++vvit) {
      Vertex_handle_pair vhp = make_vertex_handle_pair( *vvit, v);
      Classification_type status=edge_status_map[vhp];
      if (status!=INTERIOR) is_interior=false;
      if (!is_interior && status!=EXTERIOR){
        is_regular=true;
        break;
      }
  }
  
  Classification_type status=EXTERIOR;
  if (is_interior)
    status=v->is_on_chull() ? REGULAR : INTERIOR;
  else{
    if (is_regular)
      status=REGULAR;
    else{
      if ( is_vertex_Gabriel(v,Weighted_tag()) ){
        if ( is_gabriel_simplex_in_alpha_complex(v) )        status=SINGULAR;
      }
    }
  }
  v->set_classification_type(status);
}
  
template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::update_vertex_status(Vertex_handle v,Classification_type edge_status){
  Classification_type status=v->get_classification_type();
  switch(status){
    case SINGULAR:
      switch(edge_status){
        case INTERIOR: 
          status=INTERIOR;
        break;
        case EXTERIOR:
          status=EXTERIOR;
        break;
        case REGULAR:
        case SINGULAR:
          status=REGULAR;
      }
    break;
    case INTERIOR:
      switch(edge_status){
        case INTERIOR: 
        break;
        case EXTERIOR:
        case REGULAR:
        case SINGULAR:
          status=REGULAR;
      }
    break;      
    case EXTERIOR:      
      switch(edge_status){
        case EXTERIOR:
        break;
        case INTERIOR:
        case REGULAR:
        case SINGULAR:
          status=REGULAR;
      }
    break;
    case REGULAR:  
    break;
  }
  v->set_classification_type(status);
}

template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::finalize_status_of_vertex(Vertex_handle v)
{
  Classification_type status=v->get_classification_type();
  if (v->is_on_chull() && status==INTERIOR){
    v->set_classification_type(REGULAR);
    return;
  }
  if (status==INTERIOR || status==REGULAR)
    return;
  
  //when dimension is 3 any vertex has at least one incident edge,
  // thus can't be SINGULAR again (because of update_vertex_status behavior)
  CGAL_assertion(v->get_classification_type()==EXTERIOR);
  
  if ( is_vertex_Gabriel(v,Weighted_tag()) && is_gabriel_simplex_in_alpha_complex(v) )
    v->set_classification_type(SINGULAR);
}

template <class Dt>
void 
Fixed_alpha_shape_3<Dt>::initialize_status_of_vertices()
{
  #if 1 //approach avoiding extensive use of the map on 3hfli we move from 0.110983 to 0.082987
  for( Finite_vertices_iterator vit = finite_vertices_begin(); vit != finite_vertices_end();	 ++vit) 
    vit->set_classification_type(SINGULAR);
  for (typename Edge_status_map::const_iterator eit=edge_status_map.begin();eit!=edge_status_map.end();++eit){
    Classification_type edge_status=eit->second;
    Vertex_handle_pair vhp=eit->first;
    update_vertex_status(vhp.first,edge_status);
    update_vertex_status(vhp.second,edge_status);
  }
  for( Finite_vertices_iterator vit = finite_vertices_begin(); vit != finite_vertices_end();	 ++vit)
    finalize_status_of_vertex(vit);
  #else
  //This method is slower because it always makes queries in the edge classification map
  for( Finite_vertices_iterator vit = finite_vertices_begin(); 
       vit != finite_vertices_end();	 ++vit) 
    set_vertex_status(vit);
  #endif
}

//---------------------------------------------------------------------

template <class Dt>
std::ostream& operator<<(std::ostream& os,  const Fixed_alpha_shape_3<Dt>& A)
  // Inserts the alpha shape into the stream `os' as an indexed face set. 
  // Precondition: The insert operator must be defined for `Point'
{
  typedef Fixed_alpha_shape_3<Dt>                  AS;
  typedef typename AS::Vertex_handle         Vertex_handle;
  typedef typename AS::Cell_handle           Cell_handle;
  typedef typename AS::Alpha_shape_vertices_iterator 
                                             Alpha_shape_vertices_iterator;
  typedef typename AS::Alpha_shape_facets_iterator
                                             Alpha_shape_facets_iterator;

  Unique_hash_map< Vertex_handle, int > V;
  int number_of_vertices = 0;

  Alpha_shape_vertices_iterator vit;
  for( vit = A.alpha_shape_vertices_begin();
       vit != A.alpha_shape_vertices_end();
       ++vit) {
    V[*vit] = number_of_vertices++;
    os << (*vit)->point() << std::endl;
  }

  Cell_handle c;
  int i;
  Alpha_shape_facets_iterator fit;
  for( fit = A.alpha_shape_facets_begin();
       fit != A.alpha_shape_facets_end();
       ++fit) {
    c = fit->first;
    i = fit->second;
    // the following ensures that regular facets are output
    // in ccw order
    if (A.classify(*fit) == AS::REGULAR && (A.classify(c) == AS::INTERIOR)){
      c = c->neighbor(i);
      i = c->index(fit->first);
    }
    int i0 = Triangulation_utils_3::vertex_triple_index(i,0);
    int i1 = Triangulation_utils_3::vertex_triple_index(i,1);
    int i2 = Triangulation_utils_3::vertex_triple_index(i,2);
    os << V[c->vertex(i0)] << ' ' 
       << V[c->vertex(i1)] << ' ' 
       << V[c->vertex(i2)] << std::endl;
  }
  return os;
}

} //namespace CGAL

#endif //CGAL_FIXED_ALPHA_SHAPE_3_H
