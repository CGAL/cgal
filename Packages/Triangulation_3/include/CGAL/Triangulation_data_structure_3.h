// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_data_structure.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================
//
// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_H

#include <pair.h>
#include <CGAL/triple.h>
#include <list.h>
#include <vector.h>
#include <map.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names.h>

#include <CGAL/Triangulation_base_vertex.h>
#include <CGAL/Triangulation_base_cell.h>

#include <CGAL/Triangulation_ds_cell.h>
#include <CGAL/Triangulation_ds_vertex.h>

#include <CGAL/Triangulation_ds_iterators.h>
#include <CGAL/Triangulation_ds_circulators.h>

template <class Tds>
class CGAL_Triangulation_ds_cell_iterator;
template <class Tds>
class CGAL_Triangulation_ds_facet_iterator;
template <class Tds>
class CGAL_Triangulation_ds_vertex_iterator;
template <class Tds>
class CGAL_Triangulation_ds_cell_circulator;

template <class Vb, class Cb>
class CGAL_Triangulation_data_structure
{

  friend istream& operator>> CGAL_NULL_TMPL_ARGS
  (istream&, CGAL_Triangulation_data_structure<Vb,Cb>&);

  friend void CGAL_Triangulation_ds_cell<Vb,Cb>::add_list
  (CGAL_Triangulation_data_structure<Vb,Cb>&);

  friend class CGAL_Triangulation_ds_cell_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_facet_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_edge_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_vertex_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  
public:

  typedef CGAL_Triangulation_ds_vertex<Vb,Cb> Vertex;
  typedef CGAL_Triangulation_ds_cell<Vb,Cb> Cell;
  typedef pair<Cell*, int>  Facet;
  typedef CGAL_triple<Cell*, int, int> Edge;

  typedef CGAL_Triangulation_data_structure<Vb,Cb> Tds;

  typedef CGAL_Triangulation_ds_cell_iterator<Tds> Cell_iterator;
  typedef CGAL_Triangulation_ds_facet_iterator<Tds> Facet_iterator;
  typedef CGAL_Triangulation_ds_edge_iterator<Tds> Edge_iterator;
  typedef CGAL_Triangulation_ds_vertex_iterator<Tds> Vertex_iterator;
  typedef CGAL_Triangulation_ds_cell_circulator<Tds> Cell_circulator;

  // CONSTRUCTORS

  inline
  CGAL_Triangulation_data_structure() 
    : _dimension(-2), _number_of_vertices(0), _list_of_cells() 
  {}

  // advanced - TBD
//   CGAL_Triangulation_data_structure(Vertex * v)
//     : _list_of_cells()
//   {
//     init(v);
//     CGAL_triangulation_postcondition( is_valid() );
//   }

  inline
  CGAL_Triangulation_data_structure(const Tds & tds)
    : _number_of_vertices(0), _list_of_cells()
    // _number_of_vertices is set to 0 so that clear() in copy_tds() works
  {
    copy_tds(tds);
  }

  // DESTRUCTOR

  ~CGAL_Triangulation_data_structure()
  {
    clear();
  }
  
  // ASSIGNEMENT

  inline
  Tds & operator= (const Tds & tds)
  {
    copy_tds(tds);
    return *this;
  }  
    
  // ACCESS FUNCTIONS

  inline 
  int number_of_vertices() const {return _number_of_vertices;}
  
  inline 
  int dimension() const {return _dimension;}

  // USEFUL CONSTANT TIME FUNCTIONS

  //  int number_of_cells() const {  }

  // SETTING
  // to be protected ?

  inline 
  void set_number_of_vertices(int n) { _number_of_vertices = n; }

  inline 
  void set_dimension(int n) { _dimension = n; }

  // MODIFY

//    void flip(Cell* f, int i)
//     {
//    }

  //INSERTION

  void insert_outside_affine_hull(Vertex* v, // new vertex
				  Vertex* star = NULL,
				  bool reorient = false) 
    // star = vertex from which we triangulate the facet of the incremented dimension
    // ( geometrically : star = infinite vertex )
    // = Null only used to insert the 1st vertex (dimension -2 to dimension -1)

    // changes the dimension

    // if (reorient) the orientation of the cells is modified
  {
    CGAL_triangulation_precondition( v != NULL );

    Cell* c;
    Cell* d;
    Cell* e;
    int i, j;

    set_number_of_vertices( number_of_vertices()+1 );
    set_dimension( dimension()+1 );
    // this is set before the switch, so that it becomes allowed to reorient
    // new facets or cells by iterating on them (otherwise the dimension is to small

    switch ( dimension() ) {

    case -1:
      // insertion of the first vertex
      // ( geometrically : infinite vertex )
      {
	c = new Cell( *this,
		      v, NULL, NULL, NULL,
		      NULL, NULL, NULL, NULL );
	v->set_cell(c);
	break;
      }

    case 0:
      // insertion of the second vertex
      // ( geometrically : first finite vertex )
      {
	CGAL_triangulation_precondition( star != NULL );
	d = new Cell( *this,
		      v, NULL, NULL, NULL,
		      star->cell(), NULL, NULL, NULL );
	v->set_cell(d);
	star->cell()->set_neighbor(0,d);
	break;
      }

    case 1:
      // insertion of the third vertex
      // ( geometrically : second finite vertex )
      {
	CGAL_triangulation_precondition( star != NULL );
	c = star->cell();
	d = c->neighbor(0);
	// the following code could be shortened :
	// if (reorient) { i=0; j=1 }
	// else { i=1; j=0 }
	// and then use i and j instead of 0 and 1
	if (reorient) {
	  c->set_vertex(0,d->vertex(0));
	  c->set_vertex(1,star);
	  c->set_neighbor(1,d);
	  d->set_vertex(1,d->vertex(0));
	  d->set_vertex(0,v);
	  d->set_neighbor(0,c);
	  e = new Cell( *this,
			star, v, NULL, NULL,
			d, c, NULL, NULL );
	  c->set_neighbor(0,e);
	  d->set_neighbor(1,e);
	}
	else {
	  c->set_vertex(1,d->vertex(0));
	  d->set_vertex(1,v);
	  d->set_neighbor(1,c);
	  e = new Cell( *this,
			v, star, NULL, NULL,
			c, d, NULL, NULL );
	  c->set_neighbor(1,e);
	  d->set_neighbor(0,e);
	}
	
	v->set_cell(d);
	break;
      }

    case 2:
      // general case : 4th vertex ( geometrically : 3rd finite vertex )
      // degenerate cases geometrically : 1st non collinear vertex
      {
	CGAL_triangulation_precondition( star != NULL );
	
	c = star->cell();
	i = c->index(star); // i== 0 or 1
	j = (1-i);
	d = c->neighbor(j);
	
	c->set_vertex(2,v);

	e = c->neighbor(i);
	Cell* cnew = c;
	Cell* enew;
	
	while( e != d ){
	  enew = new Cell( *this );
	  enew->set_vertex(i,e->vertex(j));
	  enew->set_vertex(j,e->vertex(i));
	  enew->set_vertex(2,star);
	  
	  enew->set_neighbor(i,cnew);
	  cnew->set_neighbor(j,enew); 
	  // false at the first iteration of the loop where it should be neighbor 2
	  // it is corrected after the loop
	  enew->set_neighbor(2,e);
	  // neighbor j will be set during next iteration of the loop
	  
	  e->set_vertex(2,v);
	  e->set_neighbor(2,enew);
	  
	  c = e;
	  e = e->neighbor(i);
	  cnew = enew;
	}
	
	d->set_vertex(2,v);
	d->set_neighbor(2,enew);
	enew->set_neighbor(j,d);
	
	// corrections for star->cell() :
	c = star->cell();
	c->set_neighbor(2,c->neighbor(i)->neighbor(2));
	c->set_neighbor(j,d);
	
	v->set_cell(d);
	
	if (reorient) {
	  // reorientation of all the cells
	  Vertex* vtmp;
	  Cell* ctmp;
	  Facet_iterator fit = facets_begin();
	  
	  while(fit != facets_end()) {
	    vtmp = (*fit).first->vertex(1);
	    (*fit).first->set_vertex(1,(*fit).first->vertex(0));
	    (*fit).first->set_vertex(0,vtmp);
	    
	    ctmp = (*fit).first->neighbor(1);
	    (*fit).first->set_neighbor(1,(*fit).first->neighbor(0));
	    (*fit).first->set_neighbor(0,ctmp);
	    
	    ++fit;
	  }
	}
	break;
      }

    case 3:
      // general case : 5th vertex ( geometrically : 4th finite vertex )
      // degenerate cases : geometrically 1st non coplanar vertex
      {
	CGAL_triangulation_precondition( star != NULL );

	Cell* old_cells = list_of_cells()._next_cell; 
	// used to store the beginning of the list of cells,
	// which will be past end for the list of new cell
	// in order to be able to traverse only the new cells 
	// to find the missing neighbors (we know that new Cell() puts
	// each new cell at the beginning of the list).
	
	Cell* cnew;
	Cell_iterator it = cells_begin(); 
	// allowed since the dimension has already been set to 3

	v->set_cell(&(*it)); // ok since there is at list one ``cell''
	while (it != cells_end()) {
	  it->set_vertex(3,v);
	  if ( ! it->has_vertex(star) ) {
	    cnew = new Cell( *this,
			     it->vertex(0),it->vertex(2),
			     it->vertex(1),star,
			     NULL,NULL,NULL,&(*it));
	    it->set_neighbor(3,cnew);
	  }
	  ++it;
	}

	it = cells_begin(); 
	Cell* n;
	Cell* c;
	// traversal of the new cells only, to add missing neighbors
	while ( &(*it) != old_cells ) {
	  n = it->neighbor(3); // opposite to star
	  for ( int i=0; i<3; i++ ) {
	    int j;
	    if ( i==0 ) j=0;
	    else j=3-i; // vertex 1 and vertex 2 are always switched when
	                // creating a new cell (see above)
	    if ( ( c = n->neighbor(i)->neighbor(3) ) != NULL ) {
	      // i.e. star is not a vertex of n->neighbor(i)
	      it->set_neighbor(j,c);
	      // opposite relation will be set when it arrives on c
	      // this avoids to look for the correct index 
	      // and to test whether *it already has neighbor i
	    }
	    else {
	      // star is a vertex of n->neighbor(i)
	      it->set_neighbor(j,n->neighbor(i));
	      n->neighbor(i)->set_neighbor(3,&(*it)); // neighbor opposite to v
	    }
	  }
	  ++it;
	}
	
	// reorientation of all the cells
	if (reorient) {
	  Vertex* vtmp;
	  Cell* ctmp;
	  it = cells_begin();
	  
	  while ( it != cells_end() ) {
	    vtmp = it->vertex(1);
	    it->set_vertex(1,it->vertex(0));
	    it->set_vertex(0,vtmp);

	    ctmp = it->neighbor(1);
	    it->set_neighbor(1,it->neighbor(0));
	    it->set_neighbor(0,ctmp);
	    
	    ++it;
	  }
	}
      }
    }// end switch
    
  } 
  // end insert_outside_affine_hull

  void insert_in_edge(Vertex* v, Cell* c, int i, int j)
    // inserts v in the edge of cell c with vertices i and j
  {
    CGAL_triangulation_precondition( v != NULL && c != NULL ); 
    CGAL_triangulation_precondition( i != j );
    CGAL_triangulation_precondition( dimension() >= 1 );

    Cell* cnew;
    Cell* dnew;

    switch ( dimension() ) {

    case 1:
      {
	CGAL_triangulation_precondition( (i==0 || i==1) && (j==0 || j==1) );
	cnew = new Cell(*this,
			v,c->vertex(1),NULL,NULL,
			c->neighbor(0),c,NULL,NULL);
	c->vertex(1)->set_cell(cnew);
	c->set_vertex(1,v);
	c->neighbor(0)->set_neighbor(1,cnew);
	c->set_neighbor(0,cnew);

	v->set_cell(cnew); 
	break;
      }

    case 2:
      {
	CGAL_triangulation_precondition( i>=0 && i<=2 && j>=0 && j<=2 );
	int k=3-i-j; // index of the third vertex of the facet
	Cell* d = c->neighbor(k);
	int kd = d->index(c);
	int id = d->index(c->vertex(i));
	int jd = d->index(c->vertex(j));

	cnew = new Cell(*this);
	cnew->set_vertex(i,c->vertex(i)); 
	c->vertex(i)->set_cell(cnew);
	cnew->set_vertex(j,v);
	cnew->set_vertex(k,c->vertex(k));
	c->set_vertex(i,v);

	dnew = new Cell(*this);
	dnew->set_vertex(id,d->vertex(id));
	// d->vertex(id)->cell() is cnew OK
	dnew->set_vertex(jd,v);
	dnew->set_vertex(kd,d->vertex(kd));
	d->set_vertex(id,v);

	cnew->set_neighbor(i,c);
	Cell* nj = c->neighbor(j);
	cnew->set_neighbor(j,nj);
	nj->set_neighbor(nj->index(c),cnew);
	c->set_neighbor(j,cnew);
	cnew->set_neighbor(k,dnew);

	dnew->set_neighbor(id,d);
	nj = d->neighbor(jd);
	dnew->set_neighbor(jd,nj);
	nj->set_neighbor(nj->index(d),dnew);
	d->set_neighbor(jd,dnew);
	dnew->set_neighbor(kd,cnew);

	v->set_cell(cnew);
	break;
      }

    case 3:
      {
	CGAL_triangulation_precondition( i>=0 && i<=3 && j>=0 && j<=3 );
	Vertex* vi=c->vertex(i);
	Vertex* vj=c->vertex(j);
	
	cnew = new Cell(*this, c);
	c->set_vertex(j,v);
	vj->set_cell(cnew);
	v->set_cell(c);
	c->neighbor(i)->set_neighbor(c->neighbor(i)->index(c),cnew);
	c->set_neighbor(i,cnew);
	cnew->set_vertex(i,v);
	cnew->set_neighbor(j,c);

	// the code here duplicates a large part of the code 
	// of CGAL_Triangulation_ds_cell_circulator

	int k=Cell_circulator::other(i,j);

	Cell* ctmp = c->neighbor(k);
	Cell* cprev = c;
	Cell* cnewprev = cnew;

	while ( ctmp != c ) {
	  // the current cell is duplicated. vertices and neighbors i and j
	  // are updated during the traversal.
	  // uses the field prev of the circulator
	  i = ctmp->index(vi);
	  j = ctmp->index(vj);
	  cnew = new Cell(*this, ctmp);
	  // v will become vertex j of c
	  // and vertex i of cnew
	  ctmp->set_vertex(j,v);
	  ctmp->neighbor(i)->set_neighbor(ctmp->neighbor(i)->index(ctmp),cnew);
	  ctmp->set_neighbor(i,cnew);
	  cnew->set_vertex(i,v);
	  cnew->set_neighbor(j,ctmp);

	  // neighbor relations of all cells are used
	  // to find relations between new cells
	  cnew->set_neighbor(ctmp->index(cprev),cnewprev);
	  cnewprev->set_neighbor(cprev->index(ctmp),cnew);

	  cnewprev = cnew;
	  k=Cell_circulator::other(i,j);
	  if ( ctmp->neighbor(k) == cprev ) {
	    cprev = ctmp;
	    ctmp = ctmp->neighbor(6-i-j-k);
	  }
	  else {
	    cprev = ctmp;
	    ctmp = ctmp->neighbor(k);
	  }
	}
	cnew = c->neighbor(c->index(vi));
	cnew->set_neighbor(c->index(cprev),cnewprev);
	cnewprev->set_neighbor(cprev->index(c),cnew);
	break;
      }
    }
    set_number_of_vertices(number_of_vertices() +1);
  }// end insert_in_edge

  void insert_in_facet(Vertex* v, Cell* c, int i)
    // inserts v in the facet opposite to vertex i of cell c
  {
    CGAL_triangulation_precondition( (v != NULL) && (c != NULL)); 
    CGAL_triangulation_precondition( dimension() >= 2 );

    switch ( dimension() ) {

    case 2:
      {
	CGAL_triangulation_precondition( i == 3 );
	Cell* n = c->neighbor(2);
	Cell* cnew = new Cell(*this,
			      c->vertex(0),c->vertex(1),v,NULL,
			      c, NULL,n,NULL);
	n->set_neighbor(n->index(c),cnew);
	c->set_neighbor(2,cnew);
	c->vertex(0)->set_cell(cnew);

	n = c->neighbor(1);
	Cell* dnew = new Cell(*this,
			      c->vertex(0),v,c->vertex(2),NULL,
			      c,n,cnew,NULL);
	n->set_neighbor(n->index(c),dnew);
	c->set_neighbor(1,dnew);
	cnew->set_neighbor(1,dnew);

	c->set_vertex(0,v);
	v->set_cell(c);
	break;
      }
    case 3:
      {
	CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2 || i == 3 );
	// c will be modified to have v replacing vertex(i+3)
	int i1,i2,i3;

	if ( i&1 == 0 ) {
	  i1=(i+1)&3; i2=(i+2)&3; i3=6-i-i1-i2;
	}
	else {
	  i1=(i+1)&3; i2=(i+3)&3; i3=6-i-i1-i2;
	}
	// i,i1,i2,i3 is well oriented
	// so v will "replace" the vertices in this order
	// when creating the new cells one after another from c

	Vertex* vi=c->vertex(i);
	Vertex* v1=c->vertex(i1); 
	Vertex* v2=c->vertex(i2);
	Vertex* v3=c->vertex(i3);

	// new cell with v in place of i1
	Cell* nc = c->neighbor(i1);
	Cell* cnew1 = new Cell(*this,
			       vi,v,v2,v3,
			       NULL,nc,NULL,c);
	nc->set_neighbor(nc->index(c),cnew1);
	c->set_neighbor(i1,cnew1);

	v3->set_cell(cnew1);

	// new cell with v in place of i2
	nc = c->neighbor(i2);
	Cell* cnew2 = new Cell(*this,
			       vi,v1,v,v3,
			       NULL,cnew1,nc,c);
	nc->set_neighbor(nc->index(c),cnew2);
	c->set_neighbor(i2,cnew2);
	cnew1->set_neighbor(2,cnew2); // links to previous cell

	// v replaces i3 in c
	c->set_vertex(i3,v);

	// other side of facet containing v
	Cell* d = c->neighbor(i);
	int j = d->index(c);
	int j1=d->index(v1);// triangulation supposed to be valid
	int j2=d->index(v2);
	int j3=6-j-j1-j2;
	// then the orientation of j,j1,j2,j3 depends on the parity
	// of i-j

//	if ( (j-i)%2 == 0 ) { IDIOT !!!
	  // new cell with v in place of j1
	  Cell* nd = d->neighbor(j1);
	  Cell* dnew1 = new Cell(*this,
				 d->vertex(j),v,v3,v2,
				 cnew1,nd,d,NULL);
	  nd->set_neighbor(nd->index(d),dnew1);
	  d->set_neighbor(j1,dnew1);
	  cnew1->set_neighbor(0,dnew1);
	  
	  // new cell with v in place of j2
	  nd = d->neighbor(j2);
	  Cell* dnew2 = new Cell(*this,
				 d->vertex(j),v1,v3,v,
				 cnew2,dnew1,d,nd);
	  nd->set_neighbor(nd->index(d),dnew2);
	  d->set_neighbor(j2,dnew2);
	  cnew2->set_neighbor(0,dnew2);
	  dnew1->set_neighbor(3,dnew2);
// 	}
// 	else { IDIOT !!!
// 	  // new cell with v in place of j1
// 	  Cell* nd = d->neighbor(j1);
// 	  Cell* dnew1 = new Cell(*this,
// 			  d->vertex(j),v,v2,v3,
// 			  cnew1,nd,NULL,d);
// 	  nd->set_neighbor(nd->index(d),dnew1);
// 	  d->set_neighbor(j1,dnew1);
// 	  cnew1->set_neighbor(0,dnew1);

// 	  // new cell with v in place of j2
// 	  nd = d->neighbor(j2);
// 	  Cell* dnew2 = new Cell(*this,
// 			  d->vertex(j),v1,v,v3,
// 			  cnew2,dnew1,d,nd);
// 	  nd->set_neighbor(nd->index(d),dnew2);
// 	  d->set_neighbor(j2,dnew2);
// 	  cnew2->set_neighbor(0,dnew2);
// 	  dnew1->set_neighbor(2,dnew2);
// 	}

	// v replaces i3 in d
	d->set_vertex(j3,v);
	v->set_cell(d);

	break;
      }
    }
    set_number_of_vertices(number_of_vertices() +1);
  }
  // end insert_in_facet

  void insert_in_cell(Vertex* v, Cell* c)
    //insert in cell
  {
    CGAL_triangulation_precondition( (v != NULL) && (c != NULL));
//     c->insert_in_cell(v);

    Vertex* v0 = c->vertex(0);
    Vertex* v1 = c->vertex(1);
    Vertex* v2 = c->vertex(2);
    Vertex* v3 = c->vertex(3);

    Cell* n1 = c->neighbor(1);
    Cell* n2 = c->neighbor(2);
    Cell* n3 = c->neighbor(3);

    // c will be modified to have v,v1,v2,v3 as vertices
    Cell* c3 = new Cell(*this,v0,v1,v2,v,c,NULL,NULL,n3);
    Cell* c2 = new Cell(*this,v0,v1,v,v3,c,NULL,n2,c3);
    Cell* c1 = new Cell(*this,v0,v,v2,v3,c,n1,c2,c3);

    c3->set_neighbor(1,c1);
    c3->set_neighbor(2,c2);
    c2->set_neighbor(1,c1);

    n1->set_neighbor(n1->index(c),c1);
    n2->set_neighbor(n2->index(c),c2);
    n3->set_neighbor(n3->index(c),c3);

    c->set_vertex(0,v);
    c->set_neighbor(1,c1);
    c->set_neighbor(2,c2);
    c->set_neighbor(3,c3);

    if( v0->cell() == c  ) {  v0->set_cell(c1); }
    v->set_cell(c);
    set_number_of_vertices(number_of_vertices() +1);
  }
  // end insert_in_cell

    
  // ITERATOR METHODS

  Cell_iterator cells_begin() const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_iterator(ncthis);
  }

  Cell_iterator cells_end() const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_iterator(ncthis, 1);
  }

  Facet_iterator facets_begin() const
  {
    CGAL_triangulation_precondition( dimension() >=2 );
    Tds* ncthis = (Tds*)this;
    return Facet_iterator(ncthis);
  }

  Facet_iterator facets_end() const
  {
    CGAL_triangulation_precondition( dimension() >=2 );
    Tds* ncthis = (Tds*)this;
    return Facet_iterator(ncthis,1);
  }

  Edge_iterator edges_begin() const
  {
    CGAL_triangulation_precondition( dimension() >=1 );
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis);
  }

  Edge_iterator edges_end() const
  {
    CGAL_triangulation_precondition( dimension() >=1 );
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis,1);
  }

  Vertex_iterator vertices_begin() const
  {
    CGAL_triangulation_precondition( number_of_vertices() > 0 );
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis);
  }

  Vertex_iterator vertices_end() const
  {
    CGAL_triangulation_precondition( number_of_vertices() > 0 );
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis,1);
  }

  // CIRCULATOR METHODS

  Cell_circulator incident_cells(Edge e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_circulator(ncthis,e);
  }

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const
  {
    switch ( dimension() ) {
    case 3:
      {
	int vertex_count;
	if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
	if ( number_of_vertices() != vertex_count ) {
	  if (verbose) { cerr << "false number of vertices" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}

	int edge_count;
	if ( ! count_edges(edge_count,verbose,level) ) {return false;}
	int facet_count;
	if ( ! count_facets(facet_count,verbose,level) ) {return false;}
	int cell_count;
	if ( ! count_cells(cell_count,verbose,level) ) {return false;}

	// Euler relation 
	if ( cell_count - facet_count + edge_count - vertex_count != 0 ) {
	  if (verbose) { cerr << "Euler relation unsatisfied"<< endl; }
	  CGAL_triangulation_assertion(false); return false;
	}

	break;
      }
    case 2:
      {
	int vertex_count;
	if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
	if ( number_of_vertices() != vertex_count ) {
	  if (verbose) { cerr << "false number of vertices" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}

	int edge_count;
	if ( ! count_edges(edge_count,verbose,level) ) {return false;}
	// Euler for edges
	if ( edge_count != 3 * vertex_count - 6 ) {
	  if (verbose) { cerr << "Euler relation unsatisfied - edges/vertices" << endl;}
	  CGAL_triangulation_assertion(false); return false;
	}

	int facet_count;
	if ( ! count_facets(facet_count,verbose,level) ) {return false;}
	// Euler for facets
	if ( facet_count != 2 * vertex_count - 4 ) {
	  if (verbose) { cerr << "Euler relation unsatisfied - facets/vertices" << endl;}
	  CGAL_triangulation_assertion(false); return false;
	}
	break;
      }
    case 1:
      {
	int vertex_count;
	if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
	if ( number_of_vertices() != vertex_count ) {
	  if (verbose) { cerr << "false number of vertices" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	int edge_count;
	if ( ! count_edges(edge_count,verbose,level) ) {return false;}
	// Euler for edges
	if ( edge_count != vertex_count ) {
	  if (verbose) { cerr << "false number of edges" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	break;
      }
    case 0:
      {
	if ( number_of_vertices() < 2 ) {
	  if (verbose) { cerr << "less than 2 vertices but dimension 0" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	// no break; continue
      }
    case -1:
      {
	if ( number_of_vertices() < 1 ) {
	  if (verbose)
	    cerr << "no vertex but dimension -1" << endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	// vertex count
	int vertex_count;
	if ( ! count_vertices(vertex_count,verbose,level) )
	  return false;
	if ( number_of_vertices() != vertex_count ) {
	  if (verbose)
	    cerr << "false number of vertices" << endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
      } 
    } // end switch
    
    return true;
  } // end is_valid

  //Helping functions

  void init(Vertex*  v)
    {
    }

  void copy_tds(const Tds & tds)
  {
    map< void*, void*, less<void*> > V;
    map< void*, void*, less<void*> > F;
    Vertex*  v;
    Cell* f;

    clear();

    int n = tds.number_of_vertices();
    set_number_of_vertices(n);
    set_dimension(tds.dimension());

    if(n == 0){ return ; }

    { // create the vertices

      Vertex_iterator it=tds.vertices_begin();
      while (it != tds.vertices_end()) {
	V[&(*it)] = new Vertex( it->point() );
	++it;
      }
    }

    { // create the cells
      Cell* it = tds._list_of_cells._next_cell;
      while ( it != tds.past_end_cell() ){
	F[&(*it)]=  new Cell( *this,
			      (Vertex*) V[it->vertex(0)],
			      (Vertex*) V[it->vertex(1)],
			      (Vertex*) V[it->vertex(2)],
			      (Vertex*) V[it->vertex(3)]);
	it = it->_next_cell;
      }
    }

//    only works in dimension 3
//     { // create the cells
//       Cell_iterator it = tds.cells_begin();
//       while(it != tds.cells_end()){
// 	F[&(*it)]=  new Cell( *this,
// 			      (Vertex*) V[it->vertex(0)],
// 			      (Vertex*) V[it->vertex(1)],
// 			      (Vertex*) V[it->vertex(2)],
// 			      (Vertex*) V[it->vertex(3)]);
// 	++(it);
//       }
//     }

    { // link the vertices to a cell
      Vertex_iterator it = tds.vertices_begin();
      while(it != tds.vertices_end()) {
            v = (Vertex*) V[&(*it)];
            v->set_cell( (Cell*) F[it->cell()] );
            ++it;
        }
    }

    { // hook neighbor pointers of the cells
      Cell* it = tds._list_of_cells._next_cell;
      while ( it != tds.past_end_cell() ){
	for(int j = 0; j < 4; j++){
            f = ((Cell*) F[&(*it)]);
            f->set_neighbor(j, (Cell*) F[it->neighbor(j)] );
          }
	it = it->_next_cell;
      }
    }

//     only works in dimension 3
//     { // hook neighbor pointers of the cells
//       Cell_iterator it = tds.cells_begin();
//       while(it != tds.cells_end()){
//           for(int j = 0; j < 3; j++){
//             f = ((Cell*) F[&(*it)]);
//             f->set_neighbor(j, (Cell*) F[it->neighbor(j)] );
//           }
//           ++it;
//         }
//     }

    CGAL_triangulation_postcondition( is_valid() );
  }
 
  
  void swap(Tds &tds)
  {
    list<Vertex*> fcv = _first_collinear_vertices();
    //    Vertex*  iv = infinite_vertex();
    int     nv = number_of_vertices();
    bool lin = linear();

    //the class Geom_traits is required to have a pertinent operator=
    set_first_collinear_vertices(tds.first_collinear_vertices());
    //    set_infinite_vertex(tds.infinite_vertex());
    set_number_of_vertices(tds.number_of_vertices());
    set_linear(tds.linear());

    tds.set_first_collinear_vertices(fv);
    //    tds.set_infinite_vertex(iv);
    tds.set_number_of_vertices(nv);
    tds.set_linear(lin);
  }

  void clear()
  {

    if(number_of_vertices() == 0) {
      // the list of cells must be cleared even in this case
      Cell* it=_list_of_cells._next_cell;
      while ( it != past_end_cell() ) {
	delete it;
	// uses the destructor of ds_cell, which 
	// removes the cell from the list of cells
	it=_list_of_cells._next_cell;    
      };
      // then _list_of_cells points on itself, nothing more to do
      set_dimension(-2);
      return;
    }

    list<Vertex*> Vertices;
    {// creation of a list of all vertices
      Vertex_iterator it = vertices_begin(), done = vertices_end();
      do{
	Vertices.push_front(&(*it));
      } while(++it!=done);
    }
    
    {// deletion of the cells
      // does not use the cell iterator to work in any dimension
      Cell* it=_list_of_cells._next_cell;
      while ( it != past_end_cell() ) {
	delete it;
	// uses the destructor of ds_cell, which 
	// removes the cell from the list of cells
	it=_list_of_cells._next_cell;    
      };
      // then _list_of_cells points on itself, nothing more to do
    }

    {// deletion of the vertices
      list<Vertex*>::iterator
	it=Vertices.begin(),done=Vertices.end();
      do{
	delete *it;
      } while (++it!=done);
    }

    set_number_of_vertices(0);
    set_dimension(-2);
  }


private:
  // in dimension i, number of vertices >= i+2 
  // ( the boundary of a simplex in dimension i+1 has i+2 vertices )
  int _dimension; // 
  int _number_of_vertices;
  
  // we maintain the list of cells to be able to traverse the triangulation
  // it starts with a "foo" element that will never be removed.
  // the list is circular, the foo element being used to recognize the end
  // of the list
  Cell _list_of_cells;
  
  // ACCESS FUNCTIONS

  inline
  Cell & list_of_cells() 
    {return _list_of_cells;}
  
  inline 
  Cell* past_end_cell() const 
    {
      Tds* ncthis = (Tds *)this;
      return &( ncthis->_list_of_cells );
    } 

  // used by is-valid
  bool count_vertices(int & i, bool verbose = false, int level = 0) const
    // counts AND checks the validity
  {
    i = 0;
    Vertex_iterator it = vertices_begin();
    
    while(it != vertices_end()) {
      if ( ! it->is_valid(verbose,level) ) {
	if (verbose) { cerr << "invalid vertex" << endl; }
	CGAL_triangulation_assertion(false); return false;
      }
      ++i;
      ++it;
    }
    return true;
  } 
  
  bool count_facets(int & i, bool verbose = false, int level = 0) const
    // counts but does not check
  {
    i = 0;
    Facet_iterator it = facets_begin();
    
    while(it != facets_end()) {
//       if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
// 	if (verbose) { cerr << "invalid facet" << endl;}
// 	CGAL_triangulation_assertion(false); return false;
//       }
      ++i;
      ++it;
    }
    return true;
  }

  bool count_edges(int & i, bool verbose = false, int level = 0) const
    // counts but does not check
  {
    i = 0;
    Edge_iterator it = edges_begin();
    
    while(it != edges_end()) {
//       if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
// 	if (verbose) { cerr << "invalid edge" << endl;}
// 	CGAL_triangulation_assertion(false); return false;
//       }
      ++i;
      ++it;
    }
    return true;
  }

  bool count_cells(int & i, bool verbose = false, int level = 0) const
    // counts AND checks the validity
  {
    i = 0;
    Cell_iterator it = cells_begin();
    
    while(it != cells_end()) {
      if ( ! it->is_valid(dimension(),verbose, level) ) {
	if (verbose) { cerr << "invalid cell" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      ++i;
      ++it;
    }
    return true;
  }
  
};

template < class Vb, class Cb>
istream& operator>>
(istream& is, CGAL_Triangulation_data_structure<Vb,Cb>& tds)
{

  return is;
}


template < class Vb, class Cb>
ostream& operator<<
(ostream& os, const  CGAL_Triangulation_data_structure<Vb,Cb>  &tds)
{
  
  return os;
}


#endif CGAL_TRIANGULATION_DATA_STRUCTURE_H
