#ifndef _COMPLEX_2_IN_TRIANGULATION_3_H
#define _COMPLEX_2_IN_TRIANGULATION_3_H

#include <CGAL/circulator.h>
#include <CGAL/Union_find.h>
#include <set>
#include <map>

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

  typedef std::list<Cell_handle> Cells;

  typedef typename Facets::iterator Facets_iterator;
      
  typedef typename Cells::iterator Cells_iterator;

  typedef Const_circulator_from_container<Facets> Facet_circulator;

  typedef std::map <std::pair <Vertex_handle, Vertex_handle>, 
		    std::pair<int, std::list<Facet> > >
                                                  Edge_facet_counter;

  enum Face_type{ NOT_IN_COMPLEX, ISOLATED, BOUNDARY, REGULAR, SINGULAR};


  struct Not_in_complex {
    
    bool operator()(const Facet& f) const
    { 
      return ! f.first->is_facet_on_surface(f.second) ;
    }
  };




protected:
  Triangulation_3& tri3;
  Edge_facet_counter  edge_facet_counter;

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
  
  Facet facet_with_smallest_cell_handle(Cell_handle c, int i) const {
    Cell_handle c2 = c->neighbor(i);
    if(c2 < c){
      return  std::make_pair(c2,c2->index(c));
    }
    return std::make_pair(c,i);
  }
 public:

  // Constructors

  Complex_2_in_triangulation_3 (Triangulation_3& t3):tri3(t3) {
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

  Face_type complex_subface_type (const Facet& f) const {
    return complex_subface_type (f.first, f.second);
  }

  Face_type complex_subface_type (const Cell_handle c, const int i) const {
    return (c->is_facet_on_surface(i)) ? REGULAR : NOT_IN_COMPLEX;
  }

  Face_type complex_subface_type (const Edge& e) {
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

  Face_type complex_subface_type (const Vertex_handle v) const {
    if ( v->is_visited() ) {
      if ( is_regular(v) )
      	return REGULAR;
      else return SINGULAR;
    }
    else return NOT_IN_COMPLEX;
  }


  bool is_regular(const Vertex_handle v) const {
#ifdef USE_GRAPH
    return v->is_graph_connected();
#else 
    if(v->regular_is_cached){
      return v->regular;
    } else {
      // We have to find out if there is more than one umbrella with apex v.
      // We exploit the fact that the umbrellas do not share any edge.
      // Two facets are in the same umbrella, if they share an edge.
      // We can hence use a union find data structure to compute the sets
      // of facets that build umbrellas
      // At the end we are only interested in the number of umbrellas
      Union_find<Facet> facets;
      triangulation().incident_facets(v, filter_output_iterator(std::back_inserter(facets), Not_in_complex()));
      
      typedef std::map<Vertex_handle, Union_find<Facet>::handle>  Vertex_Set_map; 
      typedef typename Vertex_Set_map::iterator Vertex_Set_map_iterator;

      Vertex_Set_map vsmap;

      for(Union_find<Facet>::iterator it = facets.begin();
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
      v->regular = (facets.number_of_sets() == 1);
      v->regular_is_cached = true;
      return v->regular;
    }
#endif
  }

  // af : added this function as calling complex_subface_type triggers update of cache
  bool is_in_complex(Vertex_handle v) const
  {
    return v->is_visited();
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

  // computes and returns the list of incident facets of v
  // in the complex
  // af: Triangluation_3 also has such a function
  // af: Question: Is it important that the smaller cell handle is taken. What does Triangulation_3


template <typename OutputIterator>
  OutputIterator incident_facets(const Vertex_handle v, OutputIterator it) const {
    std::set<Facet> soif;

    Cells loic;
    tri3.incident_cells(v, back_inserter(loic));

    for (Cells_iterator cit = loic.begin();
	 cit != loic.end();
	 ++cit) {
      Cell_handle c = (*cit);
      int i = c->index(v);

      for (int j = 0; j < 4; j++) {
	if (i != j) {
	  soif.insert(facet_with_smallest_cell_handle(c,j));
	}
      }
    }

    // keep only facets in the complex and put eveything in a list
    //  instead of a set
    for (typename std::set<Facet>::iterator fit = soif.begin();
	 fit != soif.end();
	 ++fit) {
      Facet f = *fit;
      if ( complex_subface_type(f) != NOT_IN_COMPLEX) {
	*it = f; ++it;
      }
    }

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
    v->set_visited(true);
  }



  void set_in_complex (const Facet& f) {
    set_in_complex (f.first, f.second);
  }

  void set_in_complex (const Cell_handle c, const int i) {
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    Facet f = 
      facet_with_smallest_cell_handle(c, i);

    if (tri3.dimension() == 3) {
      // if not already in the complex
      if ( complex_subface_type (c, i) == NOT_IN_COMPLEX ) {

	c->set_surface_facet(i,true);
	c2->set_surface_facet(i2,true);

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
#ifdef USE_GRAPH
	    v->add_in_graph(f, adjacent_facets(f, v));
#endif
	  }
	}
      }
    }
    else if (tri3.dimension() == 2) {
      // if not already in the complex
      if ( complex_subface_type (c, i) == NOT_IN_COMPLEX ) {

	c->set_surface_facet(i,true);

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
#ifdef USE_GRAPH
	    v->add_in_graph(f, adjacent_facets(f, v));
#endif
            // when it was singular before it is also singular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was regular
	    if((v->regular_is_cached) && (v->regular)){
	      v->regular_is_cached = false; 
	    }
	  }
	}
      }
    }
  }

  void remove_from_complex (const Vertex_handle v) {
    v->set_visited(false);
    v->regular_is_cached = false; 
  }

  void remove_from_complex (const Facet& f) {
    remove_from_complex (f.first, f.second);
  }

  void remove_from_complex (const Cell_handle c, const int i) {
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    Facet f = 
      facet_with_smallest_cell_handle(c, i);

    if (tri3.dimension() == 3) {
      // if in the complex
      if ( complex_subface_type (c, i) != NOT_IN_COMPLEX ) {

	c->set_surface_facet(i,false);
	c2->set_surface_facet(i2,false);

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
#ifdef USE_GRAPH
	    v->remove_from_graph(f);
#endif
	    // when it was regular before it is also regular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was singular
	    if((v->regular_is_cached) && (! v->regular)){
	      v->regular_is_cached = false; 
	    }
	  }
	}
      }
    }
    else if (tri3.dimension() == 2){
      // if in the complex
      if ( complex_subface_type (c, i) != NOT_IN_COMPLEX ) {	

	c->set_surface_facet(i,false);

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
#ifdef USE_GRAPH
	    v->remove_from_graph(f);
#endif

	    // when it was regular before it is also regular now, or no longer in the complex
	    // so we only have to update the regular/singular field when it was singular
	    if((v->regular_is_cached) && (! v->regular)){
	      v->regular_is_cached = false; 
	    }
	  }
	}
      }
    }
  }
}; // end Complex_2_in_triangulation_3

CGAL_END_NAMESPACE

#endif // COMPLEX_2_IN_TRIANGULATION_3_H
