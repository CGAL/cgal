// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_data_structure_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
//
// combinatorial triangulation of the boundary of a polytope
// of dimension d in dimension d+1
// for -1 <= d <= 3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_3_H

#include <CGAL/basic.h>

#include <CGAL/triple.h>
//#include <pair.h>
#include <utility>
#include <list>
#include <map>
#include <set>
#include <vector>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/Triangulation_ds_cell_3.h>
#include <CGAL/Triangulation_ds_vertex_3.h>

#include <CGAL/Triangulation_ds_iterators_3.h>
#include <CGAL/Triangulation_ds_circulators_3.h>

#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template <class TDS>
class Triangulation_ds_cell_iterator_3;
template <class TDS>
class Triangulation_ds_facet_iterator_3;
template <class TDS>
class Triangulation_ds_edge_iterator_3;
template <class TDS>
class Triangulation_ds_vertex_iterator_3;

template <class TDS>
class Triangulation_ds_cell_circulator_3;
template <class TDS>
class Triangulation_ds_facet_circulator_3;

template <class Vb, class Cb>
class Triangulation_data_structure_3
  :public Triangulation_utils_3
{

  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream&, Triangulation_data_structure_3<Vb,Cb>&);

  friend class Triangulation_ds_cell_iterator_3
  <Triangulation_data_structure_3<Vb,Cb> >;
  friend class Triangulation_ds_facet_iterator_3
  <Triangulation_data_structure_3<Vb,Cb> >;
  friend class Triangulation_ds_edge_iterator_3
  <Triangulation_data_structure_3<Vb,Cb> >;
  friend class Triangulation_ds_vertex_iterator_3
  <Triangulation_data_structure_3<Vb,Cb> >;

   friend class Triangulation_ds_cell_circulator_3
  <Triangulation_data_structure_3<Vb,Cb> >;
   friend class Triangulation_ds_facet_circulator_3
  <Triangulation_data_structure_3<Vb,Cb> >;

public:
  typedef Triangulation_ds_vertex_3<Vb,Cb> Vertex;
  typedef Triangulation_ds_cell_3<Vb,Cb> Cell;
  typedef std::pair<Cell*, int>  Facet;
  typedef triple<Cell*, int, int> Edge;

  typedef Triangulation_data_structure_3<Vb,Cb> Tds;

  typedef Triangulation_ds_cell_iterator_3<Tds> Cell_iterator;
  typedef Triangulation_ds_facet_iterator_3<Tds> Facet_iterator;
  typedef Triangulation_ds_edge_iterator_3<Tds> Edge_iterator;
  typedef Triangulation_ds_vertex_iterator_3<Tds> Vertex_iterator;

  typedef Triangulation_ds_cell_circulator_3<Tds> Cell_circulator;
  typedef Triangulation_ds_facet_circulator_3<Tds> Facet_circulator;

  // CONSTRUCTORS

  Triangulation_data_structure_3() 
    : _dimension(-2), _number_of_vertices(0), _list_of_cells() 
  {}

  Triangulation_data_structure_3(const Vertex & v)
    : _dimension(-2), _number_of_vertices(0), _list_of_cells()
  {
    insert_increase_dimension(v);
    //    CGAL_triangulation_postcondition( is_valid() );
  }

  inline
  Triangulation_data_structure_3(const Tds & tds)
    : _number_of_vertices(0), _list_of_cells()
    // _number_of_vertices is set to 0 so that clear() in copy_tds() works
  {
    copy_tds(tds);
  }

  // DESTRUCTOR

  ~Triangulation_data_structure_3()
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

  inline 
  int number_of_cells() const 
    { 
      if ( dimension() < 3 ) return 0;

      int i=0;
      Cell_iterator it = cells_begin();
      while(it != cells_end()) {
	++i;
	++it;
      }
      return i;
    }
  
  inline 
  int number_of_facets() const
    {
      if ( dimension() < 2 ) return 0;

      int i=0;
      Facet_iterator it = facets_begin();
      while(it != facets_end()) {
	++i;
	++it;
      }
      return i;
    }

  inline 
  int number_of_edges() const
    {
      if ( dimension() < 1 ) return 0;

      int i=0;
      Edge_iterator it = edges_begin();
      while(it != edges_end()) {
	++i;
	++it;
      }
      return i;
    }

  // USEFUL CONSTANT TIME FUNCTIONS


  // SETTING
  // to be protected ?

  inline 
  void set_number_of_vertices(int n) { _number_of_vertices = n; }

  inline 
  void set_dimension(int n) { _dimension = n; }

  inline 
  Cell * create_cell() 
    { 
      Cell* c = new Cell();
      add_cell( c );
      return c; 
    }

  inline 
  Cell* create_cell( Cell* c )
    {
      Cell* cnew = new Cell(c);
      add_cell( cnew );
      return cnew; 
    }
  
  inline 
  Cell* create_cell( Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3)
    {
      Cell* cnew = new Cell(v0,v1,v2,v3);
      add_cell( cnew );
      return cnew; 
    }

  inline 
  Cell* create_cell( Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3,
		     Cell* n0, Cell* n1, Cell* n2, Cell* n3)
    {
      Cell* cnew = new Cell(v0,v1,v2,v3,n0,n1,n2,n3);
      add_cell( cnew );
      return cnew; 
    }

  // not documented
  // only used by copy_tds in the TDS class
  inline 
  Cell* create_cell( Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3,
		     const Cell & old_cell )
    {
      Cell* cnew = new Cell( v0,v1,v2,v3,old_cell );
      add_cell( cnew );
      return cnew; 
    }

  inline 
  void add_cell( Cell* c ) 
    {
      CGAL_triangulation_precondition( c != NULL );
      c->_next_cell = list_of_cells()._next_cell;
      list_of_cells()._next_cell = c;
      c->_next_cell->_previous_cell = c;
      c->_previous_cell = past_end_cell();
    }

  inline 
  void delete_cell( Cell* c ) 
    { 
      CGAL_triangulation_expensive_precondition_code
	( if ( dimension() == 3 ) {
	  CGAL_triangulation_expensive_precondition( is_cell(c) ); }
	  else {
	    if ( dimension() == 2 ) {
	      CGAL_triangulation_expensive_precondition( is_facet(c,3) ); }
	    else {
	      if ( dimension() == 1 ) {
		CGAL_triangulation_expensive_precondition( is_edge(c,0,1) ); }
	      else {
		if ( dimension() == 0 ) {
		  CGAL_triangulation_expensive_precondition
		    ( is_vertex(c->vertex(0)) ); }
	      }
	    }
	  }
	  );
	  
      c->_previous_cell->_next_cell = c->_next_cell;
      c->_next_cell->_previous_cell = c->_previous_cell;
      delete( c ); 
    }
    
  // QUERIES

  bool is_vertex(Vertex* v) const;

  bool is_edge(Cell* c, int i, int j) const;
  // returns false when dimension <1 or when indices wrong

  bool is_edge(Vertex* u, Vertex* v, 
	       Cell* & c, int & i, int & j) const;
  // returns false when dimension <1

  bool is_facet(Cell* c, int i) const;
  // returns false when dimension <2 or when indices wrong

  bool is_facet(Vertex* u, Vertex* v, Vertex* w, 
		Cell* & c, int & i, int & j, int & k) const;
  // returns false when dimension <2
      
  bool is_cell(Cell* c) const;
  // returns false when dimension <3

  bool is_cell(Vertex* u, Vertex* v, Vertex* w, Vertex* t, 
	       Cell* & c, int & i, int & j, int & k, int & l) const;
  // returns false when dimension <3

  bool is_cell(Vertex* u, Vertex* v, Vertex* w, Vertex* t) const; 
  // returns false when dimension <3

  // MODIFY

  bool flip(Facet f);
  bool flip(Cell* c, int i);
  // returns false if the facet is not flippable
  // true other wise and
  // flips facet i of cell c
  // c will be replaced by one of the new cells
  void flip_flippable(Facet f);
  void flip_flippable(Cell* c, int i);
  // flips facet i of cell c
  // c will be replaced by one of the new cells
private:
  // common to flip and filp_flippable
  void flip_really(Cell* c, int i, Cell* n, int in);
public:

  bool flip(Edge e);
  bool flip(Cell* c, int i, int j);
  // returns false if the edge is not flippable
  // true otherwise and
  // flips edge i,j of cell c
  // c will be deleted
  void flip_flippable(Edge e);
  void flip_flippable(Cell* c, int i, int j);
  // flips edge i,j of cell c
  // c will be deleted
private:
  // common to flip and filp_flippable
  void flip_really(Cell* c, int i, int j,// int next, 
		   Cell* c1, Vertex* v1, int i1, int j1, int next1,
		   Cell* c2, Vertex* v2, int i2, int j2, int next2,
		   Vertex* v3);
public:

  //INSERTION

  Vertex * insert_in_cell(const Vertex & w, Cell* c);

  Vertex * insert_in_facet(const Vertex & w, const Facet & f)
    { return insert_in_facet(w,f.first,f.second); }
  
  Vertex * insert_in_facet(const Vertex & w, Cell* c, int i);
  // inserts v in the facet opposite to vertex i of cell c
  
  Vertex * insert_in_edge(const Vertex & w, const Edge & e)   
    { return insert_in_edge(w, e.first, e.second, e.third); }
  
  Vertex * insert_in_edge(const Vertex & w, Cell* c, int i, int j);   
  // inserts v in the edge of cell c with vertices i and j

  Vertex * insert_increase_dimension(const Vertex & w, // new vertex
				     Vertex* star = NULL,
				     bool reorient = false) ;
    // star = vertex from which we triangulate the facet of the
    // incremented dimension  
    // ( geometrically : star = infinite vertex )
    // = Null only used to insert the 1st vertex (dimension -2 to dimension -1)
    // -- changes the dimension
    // -- if (reorient) the orientation of the cells is modified

  // for Delaunay :
  void star_region( std::set<void*, std::less<void*> > & region, Vertex* v,
		     Cell* c, int li );
    // region is a set of connected cells
    // c belongs to region and has facet i on the boundary of region 
    // replaces the cells in region  
    // by linking v to the boundary of region 
    
private:
  Cell* create_star( std::set<void*, std::less<void*> > & region, 
		     Vertex* v,
		     Cell* c, int li );
    // creates the cells needed by star_region

public:

  // ITERATOR METHODS

  Cell_iterator cells_begin() const
  {
    //    CGAL_triangulation_precondition( dimension() == 3 );
    if ( dimension() < 3 ) return cells_end();
    Tds* ncthis = (Tds *)this;
    return Cell_iterator(ncthis);
  }

  Cell_iterator cells_end() const
  {
    //    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_iterator(ncthis, 1);
  }

  Facet_iterator facets_begin() const
  {
    if ( dimension() < 2 ) return facets_end();
    //    CGAL_triangulation_precondition( dimension() >=2 );
    Tds* ncthis = (Tds*)this;
    return Facet_iterator(ncthis);
  }

  Facet_iterator facets_end() const
  {
    //    CGAL_triangulation_precondition( dimension() >=2 );
    Tds* ncthis = (Tds*)this;
    return Facet_iterator(ncthis,1);
  }

  Edge_iterator edges_begin() const
  {
    if ( dimension() < 1 ) return edges_end();
    //    CGAL_triangulation_precondition( dimension() >=1 );
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis);
  }

  Edge_iterator edges_end() const
  {
    //    CGAL_triangulation_precondition( dimension() >=1 );
    Tds* ncthis = (Tds*)this;
    return Edge_iterator(ncthis,1);
  }

  Vertex_iterator vertices_begin() const
  {
    if ( number_of_vertices() <= 0 ) return vertices_end();
//    CGAL_triangulation_precondition( number_of_vertices() > 0 );
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis);
  }

  Vertex_iterator vertices_end() const
  {
    //    CGAL_triangulation_precondition( number_of_vertices() > 0 );
    Tds* ncthis = (Tds*)this;
    return Vertex_iterator(ncthis,1);
  }

  // CIRCULATOR METHODS

  // cells around an edge
  Cell_circulator incident_cells(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_circulator(ncthis,e);
  }
  Cell_circulator incident_cells(Cell* ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_circulator(ncthis,ce,i,j);
  }

  Cell_circulator incident_cells(const Edge & e, Cell* start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_circulator(ncthis,e,start);
  }
  Cell_circulator incident_cells(Cell* ce, int i, int j, Cell* start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Cell_circulator(ncthis,ce,i,j,start);
  }

  //facets around an edge
    Facet_circulator incident_facets(const Edge & e) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,e);
  }
  Facet_circulator incident_facets(Cell* ce, int i, int j) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,ce,i,j);
  }
  Facet_circulator incident_facets(const Edge & e, 
				   const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,e,start);
  }
  Facet_circulator incident_facets(Cell* ce, int i, int j,
				   const Facet & start) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,ce,i,j,start);
  }
  Facet_circulator incident_facets(const Edge & e, Cell* start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,e,start,f);
  }
  Facet_circulator incident_facets(Cell* ce, int i, int j, 
				   Cell* start, int f) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Tds* ncthis = (Tds *)this;
    return Facet_circulator(ncthis,ce,i,j,start,f);
  }

  // around a vertex
  void
  incident_cells(Vertex* v, std::set<Cell*, std::less<Cell*> > & cells,
		 Cell* c = NULL ) const;

  void
  incident_vertices(Vertex* v, 
		    std::set<Vertex*, std::less<Vertex*> > & vertices,
		    Cell* c = NULL ) const;

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;


  //Helping functions
  void init(Vertex*  v)
    {
    }

  Vertex* copy_tds(const Tds & tds, Vertex* vert = NULL);
    // returns the new vertex corresponding to vert in the new tds 

  void swap(Tds & tds);

  void clear();

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

  // used by is-valid :
  bool count_vertices(int & i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
  bool count_facets(int & i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_edges(int & i, bool verbose = false, int level = 0) const;
  // counts but does not check
  bool count_cells(int & i, bool verbose = false, int level = 0) const;
  // counts AND checks the validity
  
};

template < class Vb, class Cb>
std::istream& operator >>
(std::istream& is, Triangulation_data_structure_3<Vb,Cb>& tds)
  // reads :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices 
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{

  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
  typedef  Tds::Vertex  Vertex;
  typedef  Tds::Cell Cell;
  typedef  Tds::Edge Edge;
  typedef  Tds::Facet Facet;
  //  typedef typename Vb::Point Point;

  tds.clear();

  int i;
  int n, d;
  is >> d >> n;
  tds.set_dimension(d);
  tds.set_number_of_vertices(n);

  if(n == 0){
    return is;
  }

  //  Point p;
  std::map< int, Vertex*, std::less<int> > V;
  //  std::vector<Vertex*> V(n);
  
  // creation of the vertices    

  for (i=0; i < n; i++) {
    //    is >> p;
    //    V[i] = new Vertex(p);
    V[i] = new Vertex();
  }

  std::map< int, Cell*, std::less<int> > C;
  int m;
 
  read_cells(is, tds, n, V, m, C);
  CGAL_triangulation_assertion( tds.is_valid(false) );
  return is;
}


template < class Vb, class Cb>
std::ostream& operator<<
(std::ostream& os, const Triangulation_data_structure_3<Vb,Cb>  &tds)
  // writes :
  // the dimension
  // the number of vertices
  // the number of cells
  // the cells by the indices of their vertices 
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
  typedef  Tds::Vertex  Vertex;
  typedef  Tds::Cell Cell;
  typedef  Tds::Edge Edge;
  typedef  Tds::Facet Facet;
  typedef  Tds::Vertex_iterator  Vertex_iterator;
  typedef  Tds::Cell_iterator  Cell_iterator;
  typedef  Tds::Edge_iterator  Edge_iterator;
  typedef  Tds::Facet_iterator  Facet_iterator;

  std::map< void*, int, std::less<void*> > V;
  //  std::map< void*, int, less<void*> > C;

  // outputs dimension and number of vertices
  int n = tds.number_of_vertices();
  switch ( tds.dimension() ) {
  case 3:
    {
      if(is_ascii(os)){
        os << tds.dimension() << std::endl << n << std::endl;
      } else {
        os << tds.dimension() << n;
      }
      break;
    }
  case 2:
    {
      if(is_ascii(os)){
        os << tds.dimension() << std::endl << n << std::endl;
      } else {
        os << tds.dimension() << n;
      }
      break;
    }
  case 1:
    {
      if(is_ascii(os)){
        os << tds.dimension() << std::endl << n << std::endl;
      } else {
        os << tds.dimension() << n ;
      }
      break;
    }
  case 0:
    {
      if(is_ascii(os)){
	os << tds.dimension() << std::endl << n << std::endl;
      } else {
	os << tds.dimension() << n;
      }
      break;
    }
  default:
    {
      if(is_ascii(os)){
	os << tds.dimension() << std::endl << n << std::endl;
      } else {
	os << tds.dimension() << n;
      }
    }
  }

  if (n == 0){
    return os;
  }
  
  // index the vertices
  int i = 0;
  Vertex_iterator it = tds.vertices_begin();
    
  while(it != tds.vertices_end()){
    V[&(*it)] = i++;
    //    os << it->point();
    //    if(is_ascii(os)){
    //      os << std::endl;
    //    }
    ++it;
  }
  CGAL_triangulation_assertion( i == n );

  //  print_cells(os, tds, n, V);
  print_cells(os, tds, V);

  return os;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_vertex(Vertex* v) const
    {
      Vertex_iterator it = vertices_begin();
      while (it != vertices_end()) {
	if ( v == &(*it) ) return true;
	++it;
      }
      return false;
    }

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge(Vertex* u, Vertex* v, 
	Cell* & c, int & i, int & j) const
  // returns false when dimension <1
    {
      Edge_iterator it = edges_begin();
      while ( it != edges_end() ) {
	if ( ( ((*it).first)->has_vertex(u,i) )
	     && ( ((*it).first)->has_vertex(v,j) ) ) {
	  c = (*it).first;
	  return true;
	}
	++it;
      }
      return false;
    }

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_edge(Cell* c, int i, int j) const
  // returns false when dimension <1
{
  if ( i==j ) return false;
  if ( (i<0) || (j<0) ) return false;
  if ( (dimension() == 1) && ((i>1) || (j>1)) ) return false;
  if ( (dimension() == 2) && ((i>2) || (j>2)) ) return false;
  if ((i>3) || (j>3)) return false;

  Edge_iterator it = edges_begin();
  while ( it != edges_end() ) {
    if ( (*it).first == c ) return true;
    ++it;
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_facet(Vertex* u, Vertex* v, Vertex* w, 
	 Cell* & c, int & i, int & j, int & k) const
  // returns false when dimension <2
    {
      Facet_iterator it = facets_begin();
      while ( it != facets_end() ) {
	if ( ( ((*it).first)->has_vertex(u,i) )
	     && ( ((*it).first)->has_vertex(v,j) )
	     && ( ((*it).first)->has_vertex(w,k) ) ) {
	  c = (*it).first;
	  return true;
	}
	++it;
      }
      return false;
    }

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_facet(Cell* c, int i) const
  // returns false when dimension <2
{
  if (i<0) return false;
  if ( (dimension() == 2) && (i!=3) ) return false;
  if (i>3) return false;
  Facet_iterator it = facets_begin();
  while ( it != facets_end() ) {
    if ( (*it).first == c ) return true;
    ++it;
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell( Cell* c ) const
  // returns false when dimension <3
{
  if ( c == NULL ) return false;
  Cell_iterator it = cells_begin();
  while ( it != cells_end() ) {
    if ( c == &(*it) ) {
      return true;
    }
    ++it;
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex* u, Vertex* v, Vertex* w, Vertex* t, 
	Cell* & c, int & i, int & j, int & k, int & l) const
  // returns false when dimension <3
{
  Cell_iterator it = cells_begin();
  while ( it != cells_end() ) {
    if ( ( it->has_vertex(u,i) )
	 && ( it->has_vertex(v,j) )
	 && ( it->has_vertex(w,k) ) 
	 && ( it->has_vertex(t,l) ) ) {
      c = &(*it);
      return true;
    }
    ++it;
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
is_cell(Vertex* u, Vertex* v, Vertex* w, Vertex* t) const
  // returns false when dimension <3
{
  Cell_iterator it = cells_begin();
  while ( it != cells_end() ) {
    if ( ( it->has_vertex(u) ) &&
	 ( it->has_vertex(v) ) &&
	 ( it->has_vertex(w) ) &&
	 ( it->has_vertex(t) ) ) {
      return true;
    }
    ++it;
  }
  return false;
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
flip( Facet f )
{
  return flip( f.first, f.second);
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
flip( Cell* c, int i )
  // returns false if the facet is not flippable
  // true other wise and
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4) 
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell * n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  std::set< Vertex*,std::less<Vertex*> > setc;
  incident_vertices( c->vertex(i), setc );
  if ( setc.find( n->vertex(in) ) != setc.end() ) return false;

  flip_really(c,i,n,in);
  return true;
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable( Facet f )
{
  return flip_flippable( f.first, f.second );
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable( Cell* c, int i )
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4) 
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Cell * n = c->neighbor(i);
  int in = n->index(c);

  // checks that the facet is flippable,
  // ie the future edge does not already exist
  typedef std::set< Vertex*,std::less<Vertex*> > set_of_vertices;
  CGAL_triangulation_expensive_precondition_code
    ( set_of_vertices setc; );
  CGAL_triangulation_expensive_precondition_code
    ( incident_vertices( c->vertex(i), setc ); );
  CGAL_triangulation_expensive_precondition
    ( ( setc.find( n->vertex(in) ) == setc.end() ) );

  flip_really(c,i,n,in);
}

template < class Vb, class Cb>
inline void
Triangulation_data_structure_3<Vb,Cb>::
flip_really( Cell* c, int i, Cell* n, int in )
  // private - used by flip and flip_flippable
{
  int i1 = (i+1)&3;
  int i2 = (i+2)&3;
  int i3 = (i+3)&3;

  int in1 = n->index(c->vertex(i1));
  int in2 = n->index(c->vertex(i2));
  int in3 = n->index(c->vertex(i3));

  c->set_neighbor( i, n->neighbor(in3) );
  n->neighbor(in3)->set_neighbor( n->neighbor(in3)->index(n), c );
  c->set_vertex( i3, n->vertex(in) );

  n->set_neighbor( in, c->neighbor(i1) );
  c->neighbor(i1)->set_neighbor( c->neighbor(i1)->index(c), n );
  n->set_vertex( in1, c->vertex(i) );

  Cell* cnew;
  if ( (i%2) == 0 ) 
    cnew = create_cell( c->vertex(i), c->vertex(i1),
			n->vertex(in), n->vertex(in3), 
			n->neighbor(in2), n,
			c->neighbor(i2), c );
  else
    cnew = create_cell( c->vertex(i1), c->vertex(i),
			n->vertex(in), n->vertex(in3),
			n, n->neighbor(in2), 
			c->neighbor(i2), c );

  c->neighbor(i2)->set_neighbor( c->neighbor(i2)->index(c), cnew);
  n->neighbor(in2)->set_neighbor( n->neighbor(in2)->index(n), cnew);
  
  c->set_neighbor( i1, n );
  c->set_neighbor( i2, cnew );
  n->set_neighbor( in2, cnew );
  n->set_neighbor( in3, c );

  c->vertex(i1)->set_cell(cnew);
  c->vertex(i2)->set_cell(c);
  n->vertex(in3)->set_cell(n);
  // to be implemented : 2d case
  // CGAL_triangulation_precondition( (0<=i) && (i<3) );
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
flip( Edge e )
{
  return flip( e.first, e.second, e.third );
}

template < class Vb, class Cb>
bool
Triangulation_data_structure_3<Vb,Cb>::
flip( Cell* c, int i, int j )
  // flips edge i,j of cell c
{
  CGAL_triangulation_precondition( (dimension() == 3) 
				   && (0<=i) && (i<4) 
				   && (0<=j) && (j<4)
				   && ( i != j )
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  // checks that the edge is flippable ie degree 3
  int degree = 0;
  Cell_circulator ccir = incident_cells(c,i,j);
  Cell_circulator cdone = ccir;
  do {
    ++degree;
    ++ccir;
  } while ( ccir != cdone );

  if ( degree != 3 ) return false;
  
  int next = next_around_edge(i,j);
  Cell* c1 = c->neighbor( next );
  Vertex* v1 = c->vertex( next ); // will become vertex of c1
  int i1 = c1->index( c->vertex(i) );
  int j1 = c1->index( c->vertex(j) );

  int next1 = next_around_edge(i1,j1);
  Cell* c2  = c1->neighbor( next1 );
  Vertex* v2 = c1->vertex( next1 ); // will become vertex of c2
  int i2 = c2->index( c->vertex(i) );
  int j2 = c2->index( c->vertex(j) );

  int next2 = next_around_edge(i2,j2);
  Vertex* v3 = c2->vertex( next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  if ( is_cell(v1,v2,v3,c->vertex(i)) ) return false;
  if ( is_cell(v1,v2,v3,c->vertex(j)) ) return false;

  //  flip_really(c,i,j,next,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);
  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);

  return true;
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable( Edge e )
{
  return flip_flippable( e.first, e.second, e.third );
}

template < class Vb, class Cb>
void
Triangulation_data_structure_3<Vb,Cb>::
flip_flippable( Cell* c, int i, int j )
  // flips edge i,j of cell c
{
  CGAL_triangulation_precondition( (dimension() == 3) 
				   && (0<=i) && (i<4) 
				   && (0<=j) && (j<4)
				   && ( i != j )
				   && (number_of_vertices() > 6) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  // checks that the edge is flippable ie degree 3
  CGAL_triangulation_precondition_code
    ( int degree = 0; );
  CGAL_triangulation_precondition_code
    ( Cell_circulator ccir = incident_cells(c,i,j); );
  CGAL_triangulation_precondition_code
    ( Cell_circulator cdone = ccir; );
  CGAL_triangulation_precondition_code
    ( do {
      ++degree;
      ++ccir;
    } while ( ccir != cdone ); );

  CGAL_triangulation_precondition( degree == 3 );
  
  int next = next_around_edge(i,j);
  Cell* c1 = c->neighbor( next );
  Vertex* v1 = c->vertex( next ); // will become vertex of c1
  int i1 = c1->index( c->vertex(i) );
  int j1 = c1->index( c->vertex(j) );

  int next1 = next_around_edge(i1,j1);
  Cell* c2  = c1->neighbor( next1 );
  Vertex* v2 = c1->vertex( next1 ); // will become vertex of c2
  int i2 = c2->index( c->vertex(i) );
  int j2 = c2->index( c->vertex(j) );

  int next2 = next_around_edge(i2,j2);
  Vertex* v3 = c2->vertex( next2 );

  // checks that the edge is flippable,
  // is the future cells do not already exist
  CGAL_triangulation_expensive_precondition
    ( ! is_cell(v1,v2,v3,c->vertex(i)) );
  CGAL_triangulation_expensive_precondition
    ( ! is_cell(v1,v2,v3,c->vertex(j)) );

  //  flip_really(c,i,j,next,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);
  flip_really(c,i,j,c1,v1,i1,j1,next1,c2,v2,i2,j2,next2,v3);
}

template < class Vb, class Cb>
inline void
Triangulation_data_structure_3<Vb,Cb>::
flip_really( Cell* c, int i, int j,// int next, 
	     Cell* c1, Vertex* v1, int i1, int j1, int next1,
	     Cell* c2, Vertex* v2, int i2, int j2, int next2,
	     Vertex* v3 )
{
  c->vertex(i)->set_cell(c1);
  c->vertex(j)->set_cell(c2);

  c1->set_vertex( j1, v1 );
  v1->set_cell(c1);
  c2->set_vertex( i2, v2 );
  v2->set_cell(c2);

  c1->set_neighbor( next1, c2->neighbor(j2) );
  c2->neighbor(j2)->set_neighbor( c2->neighbor(j2)->index(c2), c1 );
  c2->set_neighbor( c2->index(v1), c1->neighbor(i1) );
  c1->neighbor(i1)->set_neighbor( c1->neighbor(i1)->index(c1), c2 );

  c1->set_neighbor( i1, c2 );
  c2->set_neighbor( j2, c1 );

  c1->set_neighbor( 6-i1-j1-next1 , c->neighbor(j) );
  c->neighbor(j)->set_neighbor( c->neighbor(j)->index(c), c1 );
  c2->set_neighbor( next2, c->neighbor(i) );
  c->neighbor(i)->set_neighbor( c->neighbor(i)->index(c), c2 );

  v3->set_cell( c2 );

  delete_cell( c );

}

// not documented
template < class Vb, class Cb>
std::istream& 
read_cells(std::istream& is,
	   Triangulation_data_structure_3<Vb,Cb>  &tds,
	   int n,
	   // std::vector<void*> &V(n),
	   std::map< int, 
	             typename Triangulation_data_structure_3<Vb,Cb>::Vertex*, 
	             std::less<int> > &V,
	   // std::vector<void*> &C(m)
	   int & m,
	   std::map< int, 
                     typename Triangulation_data_structure_3<Vb,Cb>::Cell*, 
	             std::less<int> > &C )
{
  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
  typedef  Tds::Vertex  Vertex;
  typedef  Tds::Cell Cell;
  typedef  Tds::Edge Edge;
  typedef  Tds::Facet Facet;
 
    // creation of the cells and neighbors
  int i;
  switch (tds.dimension()) {
  case 3:
    {
      is >> m;
      Cell* c;

      int i0, i1, i2, i3;
      for(i = 0; i < m; i++) {
	is >> i0 >> i1 >> i2 >> i3;
	c = tds.create_cell(V[i0], V[i1], V[i2], V[i3]);
	C[i] = c;
	V[i0]->set_cell(c);
	V[i1]->set_cell(c);
	V[i2]->set_cell(c);
	V[i3]->set_cell(c);
      }
      for(i = 0; i < m; i++) {
        is >> i0 >> i1 >> i2 >> i3;
        c = C[i];
        c->set_neighbor(0, C[i0]);
        c->set_neighbor(1, C[i1]);
        c->set_neighbor(2, C[i2]);
        c->set_neighbor(3, C[i3]);
      }
      break;
    }
  case 2:
    {
      is >> m;
      Cell* c;

      int i0, i1, i2;
      for(i = 0; i < m; i++) {
	is >> i0 >> i1 >> i2;
	c = tds.create_cell(V[i0], V[i1], V[i2], NULL);
	C[i] = c;
	V[i0]->set_cell(c);
	V[i1]->set_cell(c);
	V[i2]->set_cell(c);
      }
      for(i = 0; i < m; i++) {
        is >> i0 >> i1 >> i2;
	c = C[i];
        c->set_neighbor(0, C[i0]);
        c->set_neighbor(1, C[i1]);
        c->set_neighbor(2, C[i2]);
      }
      break;
    }
  case 1:
    {
      is >> m;
     Cell* c;

      int i0, i1;
      for(i = 0; i < m; i++) {
	is >> i0 >> i1;
	c = tds.create_cell(V[i0], V[i1], NULL, NULL);
	C[i] = c;
	V[i0]->set_cell(c);
	V[i1]->set_cell(c);
      }
      for(i = 0; i < m; i++) {
        is >> i0 >> i1;
	c = C[i];
        c->set_neighbor(0, C[i0]);
        c->set_neighbor(1, C[i1]);
      }
      break;
    }
  case 0:
    {
      m = 2;
      Cell* c;

      CGAL_triangulation_assertion( (n == 2) );
      for (i=0; i < 2; i++) {
	c = tds.create_cell(V[i], NULL, NULL, NULL);
	C[i] = c;
	V[i]->set_cell(c);
      }
      for (i=0; i < 2; i++) {
	c = C[i];
        c->set_neighbor(0, C[1-i]);
      }
      break;
    }
  case -1:
    {
      m = 1;
      Cell* c;
      CGAL_triangulation_assertion( (n == 1) );
      c = tds.create_cell(V[0], NULL, NULL, NULL);
      C[0] = c;
      V[0]->set_cell(c);
      break;
    }
  }
  return is;
}
// not documented
template < class Vb, class Cb>
std::ostream& 
print_cells(std::ostream& os, 
	    const Triangulation_data_structure_3<Vb,Cb>  &tds,
	    //	    int n,
	    std::map< void*, int, std::less<void*> > &V )
{
  typedef Triangulation_data_structure_3<Vb,Cb> Tds;
  typedef  Tds::Vertex  Vertex;
  typedef  Tds::Cell Cell;
  typedef  Tds::Edge Edge;
  typedef  Tds::Facet Facet;
  typedef  Tds::Vertex_iterator  Vertex_iterator;
  typedef  Tds::Cell_iterator  Cell_iterator;
  typedef  Tds::Edge_iterator  Edge_iterator;
  typedef  Tds::Facet_iterator  Facet_iterator;

  std::map< void*, int, std::less<void*> > C;

  int i = 0;
  int j;
  int m;
  
  switch ( tds.dimension() ) {
  case 3:
    {
      m = tds.number_of_cells();
      os << m;
      if(is_ascii(os)){ os << std::endl;}

      // write the cells
      Cell_iterator it = tds.cells_begin();
      while( it != tds.cells_end() ) {
	C[&(*it)] = i++;
	for(j = 0; j < 4; j++){
	  os << V[it->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==3 ) {
	      os << std::endl;
	    } else {
	      os << ' ';
	    }
	  }
	}
	++it;
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      it = tds.cells_begin();
      while ( it != tds.cells_end() ) {
	for (j = 0; j < 4; j++) {
	  os << C[&(* it->neighbor(j))];
	  if(is_ascii(os)){
	    if(j==3) {
	      os << std::endl;
	    } else {
	      os <<  ' ';
	    }
	  }
	}
	++it;
      }
      break;
    }
  case 2:
    {
      m = tds.number_of_facets();
      os << m;
      if(is_ascii(os)){ os << std::endl;}

      // write the facets
      Facet_iterator it = tds.facets_begin();
      while( it != tds.facets_end() ) {
	C[&*((*it).first)] = i++;
	for(j = 0; j < 3; j++){
	  os << V[(*it).first->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==2 ) {
	      os << std::endl;
	    } else {
	      os <<  ' ';
	    }
	  }
	}
	++it;
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      it = tds.facets_begin();
      while ( it != tds.facets_end() ) {
	for (j = 0; j < 3; j++) {
	  os << C[&*((*it).first->neighbor(j))];
	  if(is_ascii(os)){
	    if(j==2) {
	      os << std::endl;
	    } else {
	      os <<  ' ';
	    }
	  }
	}
	++it;
      }
      break;
    }
  case 1:
    {
      m = tds.number_of_edges();
      os << m;
      if(is_ascii(os)){ os << std::endl;}

      // write the edges
      Edge_iterator it = tds.edges_begin();
      while( it != tds.edges_end() ) {
	C[&*((*it).first)] = i++;
	for(j = 0; j < 2; j++){
	  os << V[(*it).first->vertex(j)];
	  if(is_ascii(os)) {
	    if ( j==1 ) {
	      os << std::endl;
	    } else {
	      os <<  ' ';
	    }
	  }
	}
	++it;
      }
      CGAL_triangulation_assertion( i == m );
      
      // write the neighbors
      it = tds.edges_begin();
      while ( it != tds.edges_end() ) {
	for (j = 0; j < 2; j++) {
	  os << C[&*((*it).first->neighbor(j))];
	  if(is_ascii(os)){
	    if(j==1) {
	      os << std::endl;
	    } else {
	      os <<  ' ';
	    }
	  }
	}
	++it;
      }
      break;
    }
//   default:
//     {
//       os << m;
//       if(is_ascii(os)){ os << std::endl;}
//       break;
//     }
  }
  return os;
}

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex*
Triangulation_data_structure_3<Vb,Cb>::
insert_in_cell(const Vertex & w, Cell* c)
{ //insert in cell
  CGAL_triangulation_precondition( dimension() == 3 );
  CGAL_triangulation_precondition( (c != NULL) );
  CGAL_triangulation_expensive_precondition( is_cell(c) );

  Vertex * v = new Vertex(w);

  //     c->insert_in_cell(v);

  Vertex* v0 = c->vertex(0);
  Vertex* v1 = c->vertex(1);
  Vertex* v2 = c->vertex(2);
  Vertex* v3 = c->vertex(3);

  Cell* n1 = c->neighbor(1);
  Cell* n2 = c->neighbor(2);
  Cell* n3 = c->neighbor(3);

  // c will be modified to have v,v1,v2,v3 as vertices
  Cell* c3 = create_cell(v0,v1,v2,v,c,NULL,NULL,n3);
  Cell* c2 = create_cell(v0,v1,v,v3,c,NULL,n2,c3);
  Cell* c1 = create_cell(v0,v,v2,v3,c,n1,c2,c3);

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

  return v;
}
// end insert_in_cell

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex*
Triangulation_data_structure_3<Vb,Cb>::
insert_in_facet(const Vertex & w, Cell* c, int i)
{ // inserts v in the facet opposite to vertex i of cell c

  CGAL_triangulation_precondition( (c != NULL)); 
  CGAL_triangulation_precondition( dimension() >= 2 );

  Vertex * v = new Vertex(w);

  switch ( dimension() ) {

  case 3:
    {
      CGAL_triangulation_expensive_precondition( is_cell(c) );
      CGAL_triangulation_precondition( i == 0 || i == 1 || 
				       i == 2 || i == 3 );
      // c will be modified to have v replacing vertex(i+3)
      int i1,i2,i3;

      if ( (i&1) == 0 ) {
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
      Cell* cnew1 = create_cell(vi,v,v2,v3,
				NULL,nc,NULL,c);
      nc->set_neighbor(nc->index(c),cnew1);
      c->set_neighbor(i1,cnew1);

      v3->set_cell(cnew1);

      // new cell with v in place of i2
      nc = c->neighbor(i2);
      Cell* cnew2 = create_cell(vi,v1,v,v3,
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

      // new cell with v in place of j1
      Cell* nd = d->neighbor(j1);
      Cell* dnew1 = create_cell(d->vertex(j),v,v3,v2,
				cnew1,nd,d,NULL);
      nd->set_neighbor(nd->index(d),dnew1);
      d->set_neighbor(j1,dnew1);
      cnew1->set_neighbor(0,dnew1);
	  
      // new cell with v in place of j2
      nd = d->neighbor(j2);
      Cell* dnew2 = create_cell(d->vertex(j),v1,v3,v,
				cnew2,dnew1,d,nd);
      nd->set_neighbor(nd->index(d),dnew2);
      d->set_neighbor(j2,dnew2);
      cnew2->set_neighbor(0,dnew2);
      dnew1->set_neighbor(3,dnew2);

      // v replaces i3 in d
      d->set_vertex(j3,v);
      v->set_cell(d);

      break;
    }
  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_facet(c,i) );
      Cell* n = c->neighbor(2);
      Cell* cnew = create_cell(c->vertex(0),c->vertex(1),v,NULL,
			       c, NULL,n,NULL);
      n->set_neighbor(n->index(c),cnew);
      c->set_neighbor(2,cnew);
      c->vertex(0)->set_cell(cnew);

      n = c->neighbor(1);
      Cell* dnew = create_cell(c->vertex(0),v,c->vertex(2),NULL,
			       c,n,cnew,NULL);
      n->set_neighbor(n->index(c),dnew);
      c->set_neighbor(1,dnew);
      cnew->set_neighbor(1,dnew);

      c->set_vertex(0,v);
      v->set_cell(c);
      break;
    }
  }
  set_number_of_vertices(number_of_vertices() +1);

  return v;
}
// end insert_in_facet

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex*
Triangulation_data_structure_3<Vb,Cb>::
insert_in_edge(const Vertex & w, Cell* c, int i, int j)   
{ // inserts v in the edge of cell c with vertices i and j
  CGAL_triangulation_precondition( c != NULL ); 
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 );

  Vertex * v = new Vertex(w);

  Cell* cnew;
  Cell* dnew;

  switch ( dimension() ) {
    
  case 3:
    {
      CGAL_triangulation_expensive_precondition( is_cell(c) );
      CGAL_triangulation_precondition( i>=0 && i<=3 && j>=0 && j<=3 );
      Vertex* vi=c->vertex(i);
      Vertex* vj=c->vertex(j);
	
      cnew = create_cell(c);
      c->set_vertex(j,v);
      vj->set_cell(cnew);
      v->set_cell(c);
      c->neighbor(i)->set_neighbor(c->neighbor(i)->index(c),cnew);
      c->set_neighbor(i,cnew);
      cnew->set_vertex(i,v);
      cnew->set_neighbor(j,c);

      // the code here duplicates a large part of the code 
      // of Triangulation_ds_cell_circulator_3

      //	int k=Cell_circulator::other(i,j);
      //	Cell* ctmp = c->neighbor(k);
      Cell* ctmp = c->neighbor( next_around_edge(i,j) );

      Cell* cprev = c;
      Cell* cnewprev = cnew;

      while ( ctmp != c ) {
	// the current cell is duplicated. vertices and neighbors i and j
	// are updated during the traversal.
	// uses the field prev of the circulator
	i = ctmp->index(vi);
	j = ctmp->index(vj);
	cnew = create_cell(ctmp);
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
	//	  k=Cell_circulator::other(i,j);
	//	  if ( ctmp->neighbor(k) == cprev ) {
	//	    cprev = ctmp;
	//	    ctmp = ctmp->neighbor(6-i-j-k);
	//	  }
	//	  else {
	//	    cprev = ctmp;
	//	    ctmp = ctmp->neighbor(k);
	//	  }
	cprev = ctmp;
	ctmp = ctmp->neighbor( next_around_edge(i,j) );
      }
      cnew = c->neighbor(c->index(vi));
      cnew->set_neighbor(c->index(cprev),cnewprev);
      cnewprev->set_neighbor(cprev->index(c),cnew);
      break;
    }

  case 2:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
      int k=3-i-j; // index of the third vertex of the facet
      Cell* d = c->neighbor(k);
      int kd = d->index(c);
      int id = d->index(c->vertex(i));
      int jd = d->index(c->vertex(j));

      cnew = create_cell();
      cnew->set_vertex(i,c->vertex(i)); 
      c->vertex(i)->set_cell(cnew);
      cnew->set_vertex(j,v);
      cnew->set_vertex(k,c->vertex(k));
      c->set_vertex(i,v);

      dnew = create_cell();
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

  case 1:
    {
      CGAL_triangulation_expensive_precondition( is_edge(c,i,j) );
      cnew = create_cell(v,c->vertex(1),NULL,NULL,
			 c->neighbor(0),c,NULL,NULL);
      c->vertex(1)->set_cell(cnew);
      c->set_vertex(1,v);
      c->neighbor(0)->set_neighbor(1,cnew);
      c->set_neighbor(0,cnew);

      v->set_cell(cnew); 
      break;
    }
  }
  set_number_of_vertices(number_of_vertices() +1);

  return v;
}// end insert_in_edge

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex*
Triangulation_data_structure_3<Vb,Cb>::
insert_increase_dimension(const Vertex & w, // new vertex
			  Vertex* star,
			  bool reorient) 
  // star = vertex from which we triangulate the facet of the
  // incremented dimension  
  // ( geometrically : star = infinite vertex )
  // = Null only used to insert the 1st vertex (dimension -2 to dimension -1)
  // changes the dimension
  // if (reorient) the orientation of the cells is modified
{  // insert()
  Vertex * v = new Vertex(w);

  Cell* c;
  Cell* d;
  Cell* e;
  int i, j;

  switch ( dimension() ) {

  case -2:
    // insertion of the first vertex
    // ( geometrically : infinite vertex )
    {
      set_number_of_vertices( number_of_vertices()+1 );
      set_dimension( dimension()+1 );

      c = create_cell( v, NULL, NULL, NULL,
		       NULL, NULL, NULL, NULL );
      v->set_cell(c);
      break;
    }

  case -1:
    // insertion of the second vertex
    // ( geometrically : first finite vertex )
    {
      CGAL_triangulation_precondition( star != NULL );
      CGAL_triangulation_precondition( is_vertex(star) ); 
      // this precondition is not expensive when there is only one vertex!

      set_number_of_vertices( number_of_vertices()+1 );
      set_dimension( dimension()+1 );

      d = create_cell( v, NULL, NULL, NULL,
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
      CGAL_triangulation_expensive_precondition( is_vertex(star) );

      set_number_of_vertices( number_of_vertices()+1 );
      set_dimension( dimension()+1 );

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
	e = create_cell( star, v, NULL, NULL,
			 d, c, NULL, NULL );
	c->set_neighbor(0,e);
	d->set_neighbor(1,e);
      }
      else {
	c->set_vertex(1,d->vertex(0));
	d->set_vertex(1,v);
	d->set_neighbor(1,c);
	e = create_cell( v, star, NULL, NULL,
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
      CGAL_triangulation_expensive_precondition( is_vertex(star) );

      set_number_of_vertices( number_of_vertices()+1 );
      set_dimension( dimension()+1 );
      // this is set now, so that it becomes allowed to reorient
      // new facets or cells by iterating on them (otherwise the
      // dimension is to small)

      c = star->cell();
      i = c->index(star); // i== 0 or 1
      j = (1-i);
      d = c->neighbor(j);
	
      c->set_vertex(2,v);

      e = c->neighbor(i);
      Cell* cnew = c;
      Cell* enew=NULL;
	
      while( e != d ){
	enew = create_cell( );
	enew->set_vertex(i,e->vertex(j));
	enew->set_vertex(j,e->vertex(i));
	enew->set_vertex(2,star);
	  
	enew->set_neighbor(i,cnew);
	cnew->set_neighbor(j,enew); 
	// false at the first iteration of the loop where it should
	// be neighbor 2 
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
    // general case : 5th vertex ( geometrically : 4th finite vertex )
    // degenerate cases : geometrically 1st non coplanar vertex
    {
      CGAL_triangulation_precondition( star != NULL );
      CGAL_triangulation_expensive_precondition( is_vertex(star) );

      set_number_of_vertices( number_of_vertices()+1 );
      set_dimension( dimension()+1 );

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
	  cnew = create_cell( it->vertex(0),it->vertex(2),
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
    
  return v;
} // end insert_increase_dimension

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
star_region( std::set<void*, std::less<void*> > & region, 
	     Vertex* v, Cell* c, int li )
  // region is a set of connected cells
  // c belongs to region and has facet i on the boundary of region 
  // replaces the cells in region  
  // by linking v to the boundary of region
{
  CGAL_triangulation_precondition( dimension() >= 2 );
  CGAL_triangulation_precondition( region.find( (void *) c )  
				   != region.end() );
  // does not check whether region is connected 
  Cell* nouv = create_star( region, v, c, li );
  v->set_cell( nouv );
  // v->set_cell( create_star( region, v, c, li ) );
  std::set<void*, std::less<void*> >::const_iterator it;
  for( it = region.begin(); it != region.end(); ++it) {
    delete( (Cell *) *it);
  }
}

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Cell*
Triangulation_data_structure_3<Vb,Cb>::
create_star( std::set<void*, std::less<void*> > & region, 
	     Vertex* v, Cell* c, int li )
  // creates the cells needed by star_region
{
  Cell* cnew;
  if ( dimension() == 3 ) {
    int i[3];
    if ( (li%2) == 1 ) {
      i[0] = (li+1)&3;
      i[1] = (li+2)&3;
      i[2] = (li+3)&3;
    }
    else {
      i[0] = (li+2)&3;
      i[1] = (li+1)&3;
      i[2] = (li+3)&3;
    }
    cnew = create_cell( c->vertex(i[0]), c->vertex(i[1]), c->vertex(i[2]), v,
			NULL, NULL, NULL, c->neighbor(li) );
    c->neighbor(li)->set_neighbor( c->neighbor(li)->index(c), cnew);

    // look for the other three neighbors of cnew
    int j1, j2;
    Cell* cur;
    Cell* n;
    for (int ii=0; ii<3; ii++) {
      cnew->vertex(ii)->set_cell(cnew);
      // indices of the vertices of cnew such that i[ii],j1,j2,li positive
      j1 = next_around_edge(i[ii],li);
      j2 = 6-i[ii]-li-j1;
      // turn around the oriented edge j1 j2
      cur = c;
      n = c->neighbor(i[ii]);
      CGAL_triangulation_assertion( next_around_edge(j1,j2)==i[ii] );//debug
      while (true) {
	j1 = n->index( cur->vertex(j1) );
	j2 = n->index( cur->vertex(j2) );
	if ( region.find( (void*) n ) == region.end() ) { 
	  //not in conflict
	  break;
	}
	CGAL_triangulation_assertion( n != c );
	cur = n;
	n = n->neighbor( next_around_edge(j1,j2) );
      }
      // now n is outside region, cur is inside
      if ( n->neighbor( next_around_edge(j2,j1) ) == cur ) {
	// neighbor relation is reciprocical, ie
	// the cell we are looking for is not yet created
	cnew->set_neighbor(ii,create_star(region,v,cur,cur->index(n)));
	continue;
      }
      // else the cell we are looking for was already created
      cnew->set_neighbor(ii,n->neighbor( next_around_edge(j2,j1) ));
    }
    return cnew;
  } // endif dimension 3
    
  // else dimension 2 
  // i1 i2 such that v,i1,i2 positive
  int i1=ccw(li);
 // traversal of the boundary of region in ccw order to create all
  // the new facets
  Cell* bound = c;
  Vertex* v1 = c->vertex(i1);
  int ind = c->neighbor(li)->index(c); // to be able to find the
  // first cell that will be created 
  Cell* cur;
  Cell* pnew = NULL;
  do {
    cur = bound;
    // turn around v2 until we reach the boundary of region
    while ( region.find( (void*) cur->neighbor(cw(i1)) ) !=
	    region.end() ) {
      // neighbor in conflict
      cur = cur->neighbor(cw(i1));
      i1 = cur->index( v1 );
    }
    // here cur has an edge on the boundary of region
    cnew = create_cell( v, v1, cur->vertex( ccw(i1) ), NULL,
			cur->neighbor(cw(i1)), NULL, pnew, NULL);
    cur->neighbor(cw(i1))->set_neighbor
      ( cur->neighbor(cw(i1))->index(cur), cnew );
    // pnew is null at the first iteration
    v1->set_cell(cnew);
    //pnew->set_neighbor( cw(pnew->index(v1)), cnew );
    if (pnew) { pnew->set_neighbor( 1, cnew );}

    bound = cur;
    i1 = ccw(i1);
    v1 = bound->vertex(i1);
    pnew = cnew;
    //} while ( ( bound != c ) || ( li != cw(i1) ) );
  } while ( v1 != c->vertex(ccw(li)) );
  // missing neighbors between the first and the last created cells
  cur = c->neighbor(li)->neighbor(ind); // first created cell
  cnew->set_neighbor( 1, cur );
  cur->set_neighbor( 2, cnew );
  return cnew;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
incident_cells(Vertex* v, std::set<Cell*, std::less<Cell*> > & cells,
	       Cell* c ) const
{
  CGAL_triangulation_precondition( v != NULL );
  CGAL_triangulation_precondition( is_vertex(v) );

  if ( dimension() < 3 ) { return; }

  if ( c == NULL ) {
    c = v->cell();
  }
  else {
    CGAL_triangulation_precondition( c->has_vertex(v) );
  }
  if ( cells.find( c ) != cells.end() ) {
    return; // c was already found
  }
  cells.insert( c );
      
  for ( int j=0; j<4; j++ ) {
    if ( j != c->index(v) ) {
      incident_cells( v, cells, c->neighbor(j) );
    }
  }
}
  
template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
incident_vertices(Vertex* v, std::set<Vertex*, std::less<Vertex*> > & vertices,
		  Cell* c ) const
{
  CGAL_triangulation_precondition( v != NULL );
  CGAL_triangulation_precondition( is_vertex(v) );
      
  if ( number_of_vertices() < 2 ) { return; }

  if ( c == NULL ) {
    c = v->cell();
  }
  else {
    CGAL_triangulation_precondition( c->has_vertex(v) );
  }

  int d = dimension();
  int j;
  int found = 0;
  for ( j=0; j <= d; j++ ) {
    if ( j != c->index(v) ) {
      if ( vertices.find( c->vertex(j) ) == vertices.end() ) {
	vertices.insert( c->vertex(j) );
      }
      else {
	found++; // c->vertex(j) was already found 
      }
    }
  }
  if ( found == 3 ) return; // c was already visited
      
  for ( j=0; j <= d; j++ ) {
    if ( j != c->index(v) ) {
      incident_vertices( v, vertices, c->neighbor(j) );
    }
  }
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
is_valid(bool verbose, int level ) const
{ // is_valid()
  switch ( dimension() ) {
  case 3:
    {
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
      if ( number_of_vertices() != vertex_count ) {
	if (verbose) { std::cerr << "false number of vertices" 
				 << std::endl; }
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
	if (verbose) { std::cerr << "Euler relation unsatisfied"
				 << std::endl; }
	CGAL_triangulation_assertion(false); return false;
      }

      break;
    }
  case 2:
    {
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
      if ( number_of_vertices() != vertex_count ) {
	if (verbose) { std::cerr << "false number of vertices" 
				 << std::endl; }
	CGAL_triangulation_assertion(false); return false;
      }

      int edge_count;
      if ( ! count_edges(edge_count,verbose,level) ) {return false;}
      // Euler for edges
      if ( edge_count != 3 * vertex_count - 6 ) {
	if (verbose) { std::cerr 
	  << "Euler relation unsatisfied - edges/vertices" << std::endl;}
	CGAL_triangulation_assertion(false); return false;
      }

      int facet_count;
      if ( ! count_facets(facet_count,verbose,level) ) {return false;}
      // Euler for facets
      if ( facet_count != 2 * vertex_count - 4 ) {
	if (verbose) { std::cerr 
	  << "Euler relation unsatisfied - facets/vertices" << std::endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      break;
    }
  case 1:
    {
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) ) {return false;}
      if ( number_of_vertices() != vertex_count ) {
	if (verbose) { std::cerr << "false number of vertices" 
				 << std::endl; }
	CGAL_triangulation_assertion(false); return false;
      }
      int edge_count;
      if ( ! count_edges(edge_count,verbose,level) ) {return false;}
      // Euler for edges
      if ( edge_count != vertex_count ) {
	if (verbose) { std::cerr << "false number of edges" << std::endl; }
	CGAL_triangulation_assertion(false); return false;
      }
      break;
    }
  case 0:
    {
      if ( number_of_vertices() < 2 ) {
	if (verbose) { std::cerr << "less than 2 vertices but dimension 0" 
				 << std::endl; }
	CGAL_triangulation_assertion(false); return false;
      }
      // no break; continue
    }
  case -1:
    {
      if ( number_of_vertices() < 1 ) {
	if (verbose)
	  std::cerr << "no vertex but dimension -1" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      // vertex count
      int vertex_count;
      if ( ! count_vertices(vertex_count,verbose,level) )
	return false;
      if ( number_of_vertices() != vertex_count ) {
	if (verbose)
	  std::cerr << "false number of vertices" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
    } 
  } // end switch
  if (verbose) { std::cerr << "valid data structure" << std::endl; }
  return true;
} // end is_valid

template <class Vb, class Cb >
Triangulation_data_structure_3<Vb,Cb>::Vertex*
Triangulation_data_structure_3<Vb,Cb>::
copy_tds(const Tds & tds, Vertex* vert )
  // returns the new vertex corresponding to vert in the new tds 
{
  if ( vert != NULL ) {
    CGAL_triangulation_precondition( tds.is_vertex(vert) );
  }

  std::map< void*, void*, std::less<void*> > V;
  std::map< void*, void*, std::less<void*> > F;
  Vertex*  v;
  Cell* f;

  clear();

  int n = tds.number_of_vertices();
  set_number_of_vertices(n);
  set_dimension(tds.dimension());

  if(n == 0){ return vert; }

  { // create the vertices

    Vertex_iterator it=tds.vertices_begin();
    while (it != tds.vertices_end()) {
      V[&(*it)] = new Vertex( *it );
      ++it;
    }
  }

  { // create the cells
    Cell* it = tds._list_of_cells._next_cell;
    while ( it != tds.past_end_cell() ){
      // 	F[&(*it)]=  new Cell( *this,
      // 			      (Vertex*) V[it->vertex(0)],
      // 			      (Vertex*) V[it->vertex(1)],
      // 			      (Vertex*) V[it->vertex(2)],
      // 			      (Vertex*) V[it->vertex(3)]);
      // modified to keep the possible additional non combinatorial
      // information 
      F[&(*it)]=  create_cell( (Vertex*) V[it->vertex(0)],
			       (Vertex*) V[it->vertex(1)],
			       (Vertex*) V[it->vertex(2)],
			       (Vertex*) V[it->vertex(3)],
			       *it);
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
  if ( vert != NULL ) 
    return (Vertex*) V[vert];
  return NULL;
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
swap(Tds & tds)
{
  // tds and *this are supposed to be valid
  int dim = dimension();
  int nb = number_of_vertices();
  Cell *l = list_of_cells()._next_cell;
  Cell* p = list_of_cells()._previous_cell;

  set_dimension(tds.dimension());
  set_number_of_vertices(tds.number_of_vertices());

  if ( tds.list_of_cells()._next_cell == &(tds.list_of_cells()) ) {
    list_of_cells()._next_cell = list_of_cells()._previous_cell 
      = &( list_of_cells() );
  }
  else {
    _list_of_cells._next_cell = tds.list_of_cells()._next_cell;
    _list_of_cells._next_cell->_previous_cell = &(_list_of_cells);
    _list_of_cells._previous_cell = tds.list_of_cells()._previous_cell;
    _list_of_cells._previous_cell->_next_cell = &(_list_of_cells);
  }

  tds._dimension = dim;
  tds._number_of_vertices = nb;

  if ( l == &(list_of_cells()) ) {
    tds._list_of_cells._next_cell = tds._list_of_cells._previous_cell 
      = &( tds._list_of_cells );
  }
  else {
    tds._list_of_cells._next_cell = l;
    tds._list_of_cells._next_cell->_previous_cell = &(tds._list_of_cells);
    tds._list_of_cells._previous_cell = p;
    tds._list_of_cells._previous_cell->_next_cell = &(tds._list_of_cells);
  }
}

template <class Vb, class Cb >
void
Triangulation_data_structure_3<Vb,Cb>::
clear()
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

  std::list<Vertex *> Vertices;
  {// creation of a list of all vertices
    Vertex_iterator it = vertices_begin(), done = vertices_end();
    do{
      Vertices.push_back(&(*it));
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
    std::list<Vertex*>::iterator
      it=Vertices.begin(),done=Vertices.end();
    while ( it != done ) {
      delete *it;
      ++it;
    }
    //       do{
    // 	delete *it;
    //       } while (++it!=done);
  }

  set_number_of_vertices(0);
  set_dimension(-2);
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
count_vertices(int & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;
  Vertex_iterator it = vertices_begin();
    
  while(it != vertices_end()) {
    if ( ! it->is_valid(verbose,level) ) {
      if (verbose) { std::cerr << "invalid vertex" << std::endl; }
      CGAL_triangulation_assertion(false); return false;
    }
    ++i;
    ++it;
  }
  return true;
} 
template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
count_facets(int & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;
  Facet_iterator it = facets_begin();
    
  while(it != facets_end()) {
    if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
      if (verbose) { std::cerr << "invalid facet" << std::endl;}
      CGAL_triangulation_assertion(false); return false;
    }
    ++i;
    ++it;
  }
  return true;
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
count_edges(int & i, bool verbose, int level) const
  // counts but does not check
{
  i = 0;
  Edge_iterator it = edges_begin();
    
  while(it != edges_end()) {
    if ( ! (*it).first->is_valid(dimension(),verbose, level) ) {
      if (verbose) { std::cerr << "invalid edge" << std::endl;}
      CGAL_triangulation_assertion(false); return false;
    }
    ++i;
    ++it;
  }
  return true;
}

template <class Vb, class Cb >
bool
Triangulation_data_structure_3<Vb,Cb>::
count_cells(int & i, bool verbose, int level) const
  // counts AND checks the validity
{
  i = 0;
  Cell_iterator it = cells_begin();
    
  while(it != cells_end()) {
    if ( ! it->is_valid(dimension(),verbose, level) ) {
      if (verbose) { std::cerr << "invalid cell" << std::endl;}
      CGAL_triangulation_assertion(false); return false;
    }
    ++i;
    ++it;
  }
  return true;
}

CGAL_END_NAMESPACE
#endif // CGAL_TRIANGULATION_DATA_STRUCTURE_3_H
