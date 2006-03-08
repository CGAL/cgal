// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_3_H

// TODO: add the iterator 
// TODO: document the output/input function of C2T3?

#include <CGAL/circulator.h>
#include <CGAL/Union_find.h>
#include <set>
#include <map>
#include <list>

CGAL_BEGIN_NAMESPACE

template < class Tr >
class Complex_2_in_triangulation_3 {

 public:

  typedef Complex_2_in_triangulation_3 < Tr > Self;

  typedef Tr Triangulation_3;

  typedef typename Triangulation_3::Vertex_handle Vertex_handle;

  typedef typename Triangulation_3::Cell_handle Cell_handle;

  typedef typename Triangulation_3::Facet Facet;

  typedef typename Triangulation_3::Edge Edge;

  typedef std::list<Facet> Facets;

  typedef std::size_t size_type;

  typedef std::list<Cell_handle> Cells;

  typedef typename Facets::iterator Facets_iterator;

  typedef typename Cells::iterator Cells_iterator;

  typedef Const_circulator_from_container<Facets> Facet_circulator;

  typedef std::map <std::pair <Vertex_handle, Vertex_handle>,
		    std::pair<int, std::list<Facet> > >
                                                  Edge_facet_counter;

  enum Face_status{ NOT_IN_COMPLEX, ISOLATED, BOUNDARY, REGULAR, SINGULAR};


  struct Not_in_complex {

    bool operator()(const Facet& f) const
    {
      //assert(f.first < f.first->neighbor(f.second));
      return ! f.first->is_facet_on_surface(f.second) ;
    }


    bool operator()(Vertex_handle v) const
    {
      //return ! v->is_visited();
      return ! v->is_in_complex();
    }

  };




protected:
  Triangulation_3& tri3;
  Edge_facet_counter  edge_facet_counter;
  size_type m_number_of_facets;

 private:
  // computes and return an ordered pair of Vertex
  std::pair<Vertex_handle, Vertex_handle>
  make_ordered_pair(const Vertex_handle vh1, const Vertex_handle vh2) const {
    if (vh1 < vh2) {
      return std::make_pair(vh1, vh2);
    }
    else {
      return std::make_pair(vh2, vh1);
    }
  }

  Facet canonical_facet(Cell_handle c, int i) const {
    Cell_handle c2 = c->neighbor(i);
    return (c2 < c) ? std::make_pair(c2,c2->index(c)) : std::make_pair(c,i);
  }
 public:

  // Constructors

  Complex_2_in_triangulation_3 (Triangulation_3& t3) 
    : tri3(t3), m_number_of_facets(0)
  {
  }

  // Access functions

  Triangulation_3& triangulation()
  {
    return tri3;
  }

  const Triangulation_3& triangulation() const
  {
    return tri3;
  }

  Face_status face_status (const Facet& f) const {
    return face_status (f.first, f.second);
  }

  Face_status face_status (const Cell_handle c, const int i) const {
    return (c->is_facet_on_surface(i)) ? REGULAR : NOT_IN_COMPLEX;
  }

  Face_status face_status (const Edge& e) {
    typename Edge_facet_counter::iterator it =
      edge_facet_counter.find(make_ordered_pair(e.first->vertex(e.second),
						e.first->vertex(e.third)));
    if (it == edge_facet_counter.end()) return NOT_IN_COMPLEX;
    switch (it->second.first){
    case 0 : return ISOLATED;
    case 1 : return BOUNDARY;
    case 2 : return REGULAR;
    default : return SINGULAR;
    }
  }

  Face_status face_status (const Vertex_handle v) const {
    if(v->validity_mark() && ! v->is_in_complex_mark()) return NOT_IN_COMPLEX;

    //test incident edges for REUGALIRITY and count BOUNDARY edges
    typename std::list<Vertex_handle> vertices;
    triangulation().incident_facets(v, std::back_inserter(vertices));
    int number_of_boundary_incident_edges = 0; //COULD BE a Bool
    for (typename std::list<Vertex_handle>::iterator vit=vertices.begin();
	 vit != vertices.end();
	 vit++ ) {
      typename Edge_facet_counter::iterator eit =
	edge_facet_counter.find(make_ordered_pair(v, *vit) );
      if (eit != edge_facet_counter.end()) {
	if ( eit->second.first == 1) ++number_of_boundary_incident_edges;
	else if  (eit->second.first != 2)   return SINGULAR;
      }
    }

    // from now on adjacent-edges are REGULAR or BOUNDARY
    int i,j;
    union_find_of_incident_facets{v,i,j);
    if ( i == 0) { return NOT_IN_COMPLEX;}
    else if ( j > 1) {return SINGULAR};
    else {// REGULAR OR BOUNDARY
      if (number_of_boundary_incident_edges == 0) return BOUNDARY;
      return REGULAR;
    }
  }


  // This function should be called only when incident edges
  // are known to be REGULAR OR BOUNDARY
  // therefore it can set the regular_or_boundary_validity_mark
  bool is_regular_or_boundary_for_vertices(const Vertex_handle v) const {
    if(v->validity_mark()){ return v->regular_or_boundary_mark();}
    int i,j;
    union_find_of_incident_facets{v,i,j);
    v->set_regular_or_boundary_validity_mark(true);
    return (j == 1);
  }


   bool is_in_complex (const Vertex_handle v) const {
    if(v->validity_mark()){ return v->is_in_complex_mark();}
    int i,j;
    union_find_of_incident_facets{v,i,j);
    return ( ! i == 0);
  }

  // extract the subset F of facets of the complex incident to v
  // set i to the number of facets in F
  // set j to the number of connected component of the adjacency graph
  //     of F
  void union_find_of_incident_facets(Vertex_handle v, int& i, int& j) {
    Union_find<Facet> facets;
    triangulation().incident_facets( v, 
				     filter_output_iterator(
					    std::back_inserter(facets), 
					    Not_in_complex()));

    typedef std::map<Vertex_handle, 
      typename Union_find<Facet>::handle>  Vertex_Set_map;
    typedef typename Vertex_Set_map::iterator Vertex_Set_map_iterator;

    Vertex_Set_map vsmap;

    for(typename Union_find<Facet>::iterator it = facets.begin();
	it != facets.end();
	++it){
      Cell_handle ch = (*it).first;
      int i = (*it).second;
      for(int j=0; j < 3; j++){
	Vertex_handle w = ch->vertex(triangulation().vertex_triple_index(i,j));
	if(w != v){
	  Vertex_Set_map_iterator vsm_it = vsmap.find(w);
	  if(vsm_it != vsmap.end()){
	    facets.unify_sets(vsm_it->second, it);
	  } else {
	    vsmap.insert(std::make_pair(w, it));
	  }
	}
      }
    }
    
    i = facets.size(); 
    j = facets.number_of_sets();
    v->set_in_complex_mark( i > 0);
    v->set_in_complex_validity_mark( true);
    v->set_regular_or_boundary_mark( j == 1);
    // do not set regular_or_boundary_validity_mark
    // because you because thes function  may have been called
    // whith SINGULAR  incident edges.
    return;
  }

  

//   // af : added this function as calling face_status  triggers update of cache
//   bool is_in_complex(Vertex_handle v) const
//   {
//     std::cerr << "Hello guys!\n";
//     return v->is_visited();
//   }

  size_type number_of_facets() const
  {
    return m_number_of_facets;
  }

  Facet_circulator incident_facets (const Edge& e) {
    // position the circulator on the first element of the facets list
    Facets& lof =
      (edge_facet_counter[make_ordered_pair(e.first->
					    vertex(e.second),
					    e.first->
					    vertex(e.third))]).second;
  Facet_circulator fcirc(&lof);
    return fcirc;
  }


  // MY TODO : turn this function into an internal function and rename it
  // because it is not conform to what the doc says.
  // The doc says that incident_facets should return a circulator
  template <typename OutputIterator>
  OutputIterator incident_facets(const Vertex_handle v, OutputIterator it) const
  {
    // We assume that for the generated facets the Cell_handle is smaller than the opposite one
    triangulation().incident_facets(v, filter_output_iterator(it, Not_in_complex()));
    return it;
  }

  // computes and returns the list of adjacent facets of f
  // with the common Vertex_handle v
  Facets adjacent_facets (const Facet& f, const Vertex_handle v) {

    Cell_handle c = f.first;
    int i = f.second;
    int iv = c->index(v);
    Edge e[2];

    // search for the two other vertices than v in f
    int k = 0;
    for (int j = 0; j < 4; j++) {
      if ( (j != i) && (j != iv) ){
	e[k] = make_triple(c, iv, j);
	k++;
      }
    }

    Facets& lof1 =
      (edge_facet_counter[make_ordered_pair(e[0].first->
					    vertex(e[0].second),
					    e[0].first->
					    vertex(e[0].third))]).second;

    Facets& lof2 =
      (edge_facet_counter[make_ordered_pair(e[1].first->
					    vertex(e[1].second),
					    e[1].first->
					    vertex(e[1].third))]).second;

    Facets lof = typename Facets::list();

    for (Facets_iterator it = lof1.begin();
	 it != lof1.end();
	 it++) {
      lof.push_back(*it);
    }

    for (Facets_iterator it = lof2.begin();
	 it != lof2.end();
	 it++) {
      lof.push_back(*it);
    }

    assert(!lof.empty());

    lof.remove(f);

    return lof;
  }

  // Setting functions

  void set_in_complex (const Vertex_handle v) {
    v->set_in_complex(true);
  }


  void set_in_complex (const Facet& f) {
    set_in_complex (f.first, f.second);
  }

  void set_in_complex (const Cell_handle c, const int i) {
    ++m_number_of_facets;
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    Facet f = canonical_facet(c, i);

    if (tri3.dimension() == 3) {
      // if not already in the complex
      if ( face_status (c, i) == NOT_IN_COMPLEX ) {

	c->set_facet_on_surface(i,true);
	c2->set_facet_on_surface(i2,true);

	// We consider only pairs made by vertices without i
	for (int j = 0; j < 4; j++) {
	  for (int k = j + 1; k < 4; k++) {
	    if ( (i != j) && (i != k) ){
	      std::pair<Vertex_handle, Vertex_handle>
		e = make_ordered_pair(c->vertex(j),
				      c->vertex(k));
	      (edge_facet_counter[e]).first++;

	      (edge_facet_counter[e]).second.push_back(f);
	    }
	  }
	}
	// add each v of f in the complex
	// add f in graph of each of these v
	// with the appropriate connexity
	for (int j = 0; j < 4; j++) {
	  if (j != i) {
	    Vertex_handle v = c->vertex(j);
	    set_in_complex(v);
	  }
	}
      }
    }
    else if (tri3.dimension() == 2) {
      // if not already in the complex
      if ( face_status (c, i) == NOT_IN_COMPLEX ) {

	c->set_facet_on_surface(i,true);

	for (int j = 0; j < 3; j++) {
	  for (int k = j + 1; k < 3; k++) {
	    if ( (i != j) && (i != k) ){
	      std::pair<Vertex_handle, Vertex_handle>
		e = make_ordered_pair(c->vertex(j),
				      c->vertex(k));
	      (edge_facet_counter[e]).first++;

	      (edge_facet_counter[e]).second.push_back(f);
	    }
	  }
	}
	// add each v of f in the complex
	// add f in graph of each of these v
	for (int j = 0; j < 3; j++) {
	  if (j != i) {
	    Vertex_handle v = c->vertex(j);
	    set_in_complex(v);
            // when it was singular before it is also singular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was regular
	    if((v->is_status_cached()) && (v->is_regular())){
	      v->set_status_cached(false);
	    }
	  }
	}
      }
    }
  }

  void remove_from_complex (const Vertex_handle v) {
    v->set_in_complex(false);
    v->set_status_cached(false);
  }

  void remove_from_complex (const Facet& f) {
    remove_from_complex (f.first, f.second);
  }

  void remove_from_complex (const Cell_handle c, const int i) {
    --m_number_of_facets;
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    Facet f = canonical_facet(c, i);

    if (tri3.dimension() == 3) {
      // if in the complex
      if ( face_status (c, i) != NOT_IN_COMPLEX ) {

	c->set_facet_on_surface(i,false);
	c2->set_facet_on_surface(i2,false);

	// update the edge counter
	for (int j = 0; j < 4; j++) {
	  for (int k = j + 1; k < 4; k++) {
	    if ( (i != j) && (i != k) ){
	      std::pair<Vertex_handle, Vertex_handle>
		e = make_ordered_pair(c->vertex(j),
				      c->vertex(k));
	      (edge_facet_counter[e]).first--;

	      (edge_facet_counter[e]).second.remove(f);
	    }
	  }
	}
	// remove v of f in the complex
	// remove f in graph of each of these v
	for (int j = 0; j < 4; j++) {
	  if (j != i) {
	    Vertex_handle v = c->vertex(j);
	    remove_from_complex(v);
	    // when it was regular before it is also regular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was singular
	    if((v->is_status_cached()) && (! v->is_regular())){
	      v->set_status_cached(false);
	    }
	  }
	}
      }
    }
    else if (tri3.dimension() == 2){
      // if in the complex
      if ( face_status (c, i) != NOT_IN_COMPLEX ) {

	c->set_facet_on_surface(i,false);

	for (int j = 0; j < 3; j++) {
	  for (int k = j + 1; k < 3; k++) {
	    if ( (i != j) && (i != k) ){
	      std::pair<Vertex_handle, Vertex_handle>
		e = make_ordered_pair(c->vertex(j),
				      c->vertex(k));
	      (edge_facet_counter[e]).first--;

	      (edge_facet_counter[e]).second.remove(f);
	    }
	  }
	}
	////////////////////////////////////////////
	////////////////// A VERIFIER QU'IL N'Y A QUE CA !!!!!!!!!!!!
	//////////////////////////////////////////
	// remove each v of f in the complex
	// remove f in graph of each of these v
	for (int j = 0; j < 3; j++) {
	  if (j != i) {
	    Vertex_handle v = c->vertex(j);
	    remove_from_complex(v);

	    // when it was regular before it is also regular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was singular
	    if((v->is_status_cached()) && (! v->is_regular())){
	      v->set_status_cached(false);
	    }
	  }
	}
      }
    }
  }
}; // end Complex_2_in_triangulation_3

CGAL_END_NAMESPACE

#endif // CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
