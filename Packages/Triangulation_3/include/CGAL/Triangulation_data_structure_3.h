// ============================================================================
//
// $Id$
//
// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_H

#include <pair.h>
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

template <class Tds>
class CGAL_Triangulation_ds_cell_iterator;
template <class Tds>
class CGAL_Triangulation_ds_facet_iterator;
template <class Tds>
class CGAL_Triangulation_ds_vertex_iterator;

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
  friend class CGAL_Triangulation_ds_vertex_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  
public:

  typedef CGAL_Triangulation_ds_vertex<Vb,Cb> Vertex;
  typedef CGAL_Triangulation_ds_cell<Vb,Cb> Cell;
  typedef pair<Cell*, int>  Facet;

  typedef CGAL_Triangulation_data_structure<Vb,Cb> Tds;

  typedef CGAL_Triangulation_ds_cell_iterator<Tds> Cell_iterator;
  typedef CGAL_Triangulation_ds_facet_iterator<Tds> Facet_iterator;
  typedef CGAL_Triangulation_ds_vertex_iterator<Tds> Vertex_iterator;

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
    // ( geometrically : infinite vertex )
    // = Null only used to insert the 1st vertex (dimension -2 to dimension -1)

    // changes the dimension

    // if (reorient) the orientation of the cells is modified
  {
    CGAL_triangulation_precondition( v != NULL );

    Cell* c;
    Cell* d;
    Cell* e;
    int i, j;

    switch ( dimension() ) {

    case -2:
      // insertion of the first vertex
      // ( geometrically : infinite vertex )
      {
	c = new Cell( *this,
		      v, NULL, NULL, NULL,
		      NULL, NULL, NULL, NULL );
	v->set_cell(c);
	break;
      }

    case -1:
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

    case 0:
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

    case 1:
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

    case 2:
      // geometrically : insertion of the 1st non coplanar vertex
      // 5th vertex ( geometrically : 4th finite vertex )
      {
	CGAL_triangulation_precondition( star != NULL );
      }
    
    }// end switch
    
    set_number_of_vertices( number_of_vertices()+1 );
    set_dimension( dimension()+1 );
  } 
  // end insert_outside_affine_hull

  void insert_in_edge(Vertex* v, Cell* c, int i, int j)
    // inserts v in the edge of cell c with vertices i and j
  {
    CGAL_triangulation_precondition( v != NULL && c != NULL ); 
    CGAL_triangulation_precondition( i != j );
    CGAL_triangulation_precondition( dimension() >= 1 );
    if (dimension() == 1) {
      CGAL_triangulation_precondition( (i==0 || i==1) && (j==0 || j==1) );
    }
    if (dimension() == 2) {
      CGAL_triangulation_precondition( i>=0 && i<=2 && j>=0 && j<=2 );
    }
    if (dimension() == 3) {
      CGAL_triangulation_precondition( i>=0 && i<=3 && j>=0 && j<=3 );
    }

    Cell* cnew;
    Cell* dnew;

    switch ( dimension() ) {

    case 1:
      {
	cnew = new Cell(*this,
			v,c->vertex(1),NULL,NULL,
			c->neighbor(0),c,NULL,NULL);
	c->vertex(1)->set_cell(cnew);
	c->set_vertex(1,v);
	c->neighbor(0)->set_neighbor(1,cnew);
	c->set_neighbor(0,cnew);
	//      v->set_cell(cnew); plus bas
	break;
      }

    case 2:
      {
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
      
	break;
      }

    case 3:
      break;
    }
    set_number_of_vertices(number_of_vertices() +1);
    v->set_cell(cnew);
  }
  // end insert_in_edge

  void insert_in_facet(Vertex* v, Cell* c, int i)
    // inserts v in the facet opposite to vertex i of cell c
  {
    CGAL_triangulation_precondition( (v != NULL) && (c != NULL)); 
    CGAL_triangulation_precondition( dimension() >= 2 );
    if (dimension() == 2) {
      CGAL_triangulation_precondition( i == 3 );
    }   
    if (dimension() == 3) {
      CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2 || i == 3 );
    }

    Cell* cnew;
    Cell* dnew;

    switch ( dimension() ) {

    case 2:
      {
	Cell* n = c->neighbor(2);
	cnew = new Cell(*this,
			c->vertex(0),c->vertex(1),v,NULL,
			c, NULL,n,NULL);
	n->set_neighbor(n->index(c),cnew);
	c->set_neighbor(2,cnew);
	c->vertex(0)->set_cell(cnew);

	n = c->neighbor(1);
	dnew = new Cell(*this,
			c->vertex(0),v,c->vertex(2),NULL,
			c,n,cnew,NULL);
	n->set_neighbor(n->index(c),dnew);
	c->set_neighbor(1,dnew);
	cnew->set_neighbor(1,dnew);

	c->set_vertex(0,v);
	break;
      }
    case 3:
      {
	break;
      }
    }
    set_number_of_vertices(number_of_vertices() +1);
    v->set_cell(cnew);
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

    set_vertex(0,v);
    set_neighbor(1,c1);
    set_neighbor(2,c2);
    set_neighbor(3,c3);

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

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const
  {
    switch ( dimension() ) {
    case 3:
      {
      }
    case 2:
      {
	// vertex count
	int vertex_count;
	if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
	if ( number_of_vertices() != vertex_count ) {
	  if (verbose) { cerr << "false number of vertices" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	// edge count to be done
	// face count
	int facet_count;
	if ( ! count_facets(facet_count,verbose,level) ) {return false;}
	if ( facet_count != 2 * vertex_count - 4 ) {
	  if (verbose) { cerr << "Euler relation unsatisfied - facets/vertices" << endl;}
	  CGAL_triangulation_assertion(false); return false;
	}
	break;
      }
    case 1:
      // edge count to be done
      // check nb of vertices = nb of edges

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
	for(int j = 0; j < 3; j++){
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
  {
    i = 0;
    Facet_iterator it = facets_begin();
    
    while(it != facets_end()) {
      if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
	if (verbose) { cerr << "invalid facet" << endl;}
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
