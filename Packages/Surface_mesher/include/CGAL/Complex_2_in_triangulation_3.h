#ifndef _COMPLEX_2_IN_TRIANGULATION_3_H
#define _COMPLEX_2_IN_TRIANGULATION_3_H

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

 protected:
  Triangulation_3& tri3;
  Edge_facet_counter  edge_facet_counter;

 private:
  // computes and return an ordered pair of Vertex
  pair<Vertex_handle, Vertex_handle>
  make_ordered_pair(const Vertex_handle vh1, const Vertex_handle vh2) const {
    if (vh1 < vh2) {
      return std::make_pair(vh1, vh2);
    }
    else {
      return std::make_pair(vh2, vh1);
    }
  }

  // computes and return the semi-facet with the smallest cell_handle
  Facet facet_with_smallest_cell_handle(const Facet& f) const {
    Cell_handle c = f.first;
    int i = f.second;
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    
    Cell_handle cmin = c;
    int imin = i;

    if (c2 < cmin) {
      cmin = c2;
      imin = i2;
    }

    return make_pair(cmin, imin);
  }

 public:

  // Constructors

  Complex_2_in_triangulation_3 (Triangulation_3& t3):tri3(t3) {
  }

  // Access functions

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
      if ( v->is_graph_connected() )
      	return REGULAR;
      else return SINGULAR;
    }
    else return NOT_IN_COMPLEX;
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
  Facets incident_facets(const Vertex_handle v) const {
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
	  Facet f = std::make_pair(c, j);
	  soif.insert(facet_with_smallest_cell_handle(f));
	}
      }
    }

    // keep only facets in the complex and put eveything in a list
    //  instead of a set
    Facets lof = Facets();
    for (typename std::set<Facet>::iterator fit = soif.begin();
	 fit != soif.end();
	 ++fit) {
      Facet f = *fit;
      if ( complex_subface_type(f) != NOT_IN_COMPLEX) {
	lof.push_back(f);
      }
    }

    return lof;
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
      facet_with_smallest_cell_handle(std::make_pair(c, i));

    if (tri3.dimension() == 3) {
      // if not already in the complex
      if ( complex_subface_type (c, i) == NOT_IN_COMPLEX ) {

	c->set_surface_facet(i,true);
	c2->set_surface_facet(i2,true);

	// We consider only pairs made by vertices without i
	for (int j = 0; j < 4; j++) {
	  for (int k = j + 1; k < 4; k++) {
	    if ( (i != j) && (i != k) ){
	      pair<Vertex_handle, Vertex_handle> 
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
	    v->add_in_graph(f, adjacent_facets(f, v));
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
	      pair<Vertex_handle, Vertex_handle> 
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
	    v->add_in_graph(f, adjacent_facets(f, v));
	  }
	}
      }
    }
  }

  void remove_from_complex (const Vertex_handle v) {
    v->set_visited(false);
  }

  void remove_from_complex (const Facet& f) {
    remove_from_complex (f.first, f.second);
  }

  void remove_from_complex (const Cell_handle c, const int i) {
    Cell_handle c2 = c->neighbor(i);
    int i2 = c2->index(c);
    Facet f = 
      facet_with_smallest_cell_handle(std::make_pair(c, i));

    if (tri3.dimension() == 3) {
      // if in the complex
      if ( complex_subface_type (c, i) != NOT_IN_COMPLEX ) {

	c->set_surface_facet(i,false);
	c2->set_surface_facet(i2,false);

	// update the edge counter
	for (int j = 0; j < 4; j++) {
	  for (int k = j + 1; k < 4; k++) {
	    if ( (i != j) && (i != k) ){
	      pair<Vertex_handle, Vertex_handle> 
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
	    v->remove_from_graph(f);
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
	      pair<Vertex_handle, Vertex_handle> 
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
	    v->remove_from_graph(f);
	  }
	}
      }
    }
  }
}; // end Complex_2_in_triangulation_3

CGAL_END_NAMESPACE

#endif // COMPLEX_2_IN_TRIANGULATION_3_H
