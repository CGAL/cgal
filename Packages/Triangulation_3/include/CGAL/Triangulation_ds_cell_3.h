// ============================================================================
//
// $Id$
//
// cell of a triangulation of any dimension <=3
//
// ============================================================================

#ifndef CGAL_TETRAHEDRALIZATION_DS_CELL_H
#define CGAL_TETRAHEDRALIZATION_DS_CELL_H

#include <CGAL/Triangulation_short_names.h>

template < class Vb, class Cb >
class CGAL_Triangulation_ds_vertex;

template < class Vb, class Cb >
class CGAL_Triangulation_data_structure;

template < class Tds>
class CGAL_Triangulation_ds_iterator_base;
template < class Tds>
class CGAL_Triangulation_ds_cell_iterator;
template < class Tds>
class CGAL_Triangulation_ds_facet_iterator;
template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator;

template < class Vb, class Cb >
class CGAL_Triangulation_ds_cell
  : public Cb
{
  friend class CGAL_Triangulation_data_structure<Vb,Cb>;

  friend class CGAL_Triangulation_ds_iterator_base
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_cell_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_facet_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_vertex_iterator
  <CGAL_Triangulation_data_structure<Vb,Cb> >;
  

public:

  typedef CGAL_Triangulation_data_structure<Vb,Cb> Tds;
  typedef CGAL_Triangulation_ds_vertex<Vb,Cb> Vertex;
  typedef typename CGAL_Triangulation_data_structure<Vb,Cb>::Facet Facet;
  typedef CGAL_Triangulation_ds_cell<Vb,Cb> Cell;

  // CONSTRUCTORS

  // used only for initializing _list_of_cells by the constructors of
  // CGAL_Triangulation_data_structure
  // private ?
  inline
  CGAL_Triangulation_ds_cell()
    : Cb(), _previous_cell(this), _next_cell(this)
  {}

  inline
  CGAL_Triangulation_ds_cell(Tds& tds)
    : Cb()
    // builds a new cell of tds and maintains the list of cells
  { add_list(tds); }
    
  CGAL_Triangulation_ds_cell(Tds& tds,
			     Vertex* v0, Vertex* v1, 
			     Vertex* v2, Vertex* v3)
    :  Cb(v0,v1,v2,v3)
  { add_list(tds); }
    
  CGAL_Triangulation_ds_cell(Tds& tds,
			     Vertex* v0, Vertex* v1, 
			     Vertex* v2, Vertex* v3,
			     Cell* n0, Cell* n1, Cell* n2, Cell* n3)
    :  Cb(v0,v1,v2,v3,n0,n1,n2,n3)
  { add_list(tds); }


  // DESTRUCTOR

  ~CGAL_Triangulation_ds_cell()
    {
      _previous_cell->_next_cell = _next_cell;
      _next_cell->_previous_cell = _previous_cell;
      // automatically calls the destructor of the face base ?
    }

  // SETTING

  inline 
  void set_vertex(int i, Vertex* v)
  {
    Cb::set_vertex(i,v);
  }
    
  inline 
   void set_neighbor(int i, Cell* n)
  {
    Cb::set_neighbor(i,n);
  }

  //  void set_vertices() inherited
      
  inline 
  void set_vertices(Vertex* v0,
		    Vertex* v1,
		    Vertex* v2,
		    Vertex* v3)
  {
    Cb::set_vertices(v0,v1,v2,v3);
  }
    
  //   void set_neighbors() inherited
     
  inline
  void set_neighbors(Cell* n0,
		     Cell* n1,
		     Cell* n2,
		     Cell* n3)
  {
    Cb::set_neighbors(n0,n1,n2,n3);
  }

  // VERTEX ACCESS

  inline
  Vertex* vertex(int i) const
  {
    return( (Vertex*) (Cb::vertex(i)));
  } 

  inline 
  bool has_vertex(const Vertex* v) const
  {
    return (Cb::has_vertex(v));
  }
    
  inline 
  bool has_vertex(const Vertex* v, int& i) const
  {
    return (Cb::has_vertex(v,i));
  }
    
  inline 
  int index(const Vertex* v) const
  {
    return(Cb::vertex_index(v));
  }

  // NEIGHBOR ACCESS

  inline 
  Cell* neighbor(int i) const
  {
    return ((Cell*) Cb::neighbor(i));
  }
    
  inline 
  bool has_neighbor(const Cell* n) const
  {
    return (Cb::has_neighbor(n));
  }
    
  inline 
  bool has_neighbor(const Cell* n, int& i) const
  {
    return (Cb::has_neighbor(n,i));
  }
    
  inline 
  int index(const Cell* n) const
  {
    return(Cb::cell_index(n));
  }
    

  // CHECKING

  bool is_valid(int dim = 3, bool verbose = false, int level = 0) const
  {
    if ( ! Cb::is_valid() ) return false;

    switch (dim) {
      
    case -2:

    case -1:
    {
      if ( vertex(0) == NULL ) {
	if (verbose) { cerr << "vertex 0 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( vertex(1) != NULL || 
	   vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose) { cerr << "vertex 1,2 or 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( neighbor(0) != NULL || neighbor(1) != NULL ||
	   neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose) { cerr << "one neighbor != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      break;
      }

    case 0:
      {
      if ( vertex(0) == NULL ) {
	if (verbose) { cerr << "vertex 0 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( neighbor (0) == NULL ) {
	if (verbose) { cerr << "neighbor 0 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( vertex(1) != NULL || 
	   vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose) { cerr << "vertex 1, 2 or 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( neighbor(1) != NULL ||
	   neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose) { cerr << "neighbor 1, 2 or 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }

      if ( ! neighbor(0)->has_vertex(vertex(0)) ) {
	if (verbose) { cerr << "neighbor 0 does not have vertex 0" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      break;
      }

    case 1:
      {
      Vertex* v0; Vertex* v1;
      Cell* n0; Cell* n1;

      if ( v0 == NULL || v1 == NULL ) {
	if (verbose) { cerr << "vertex 0 or 1 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( n0 == NULL || n1 == NULL ) {
	if (verbose) { cerr << "neighbor 0 or 1 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose) { cerr << "vertex 2 or 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose) { cerr << "neighbor 2 or 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }

      if ( v0 !=  n0->vertex(1) ) {
	if (verbose) { cerr << 
		       "neighbor 0 does not have vertex 0 as vertex 1" 
			    << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( v1 != n1->vertex(0) ) {
	if (verbose) { cerr << 
			 "neighbor 1 does not have vertex 1 as vertex 0" 
			    << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      
      if ( this != n0->neighbor(1) ) {
	if (verbose) { cerr << 
			 "neighbor 0 does not have this as neighbor 1" 
			    << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( this != n1->neighbor(0) ) {
	if (verbose) { cerr << 
			 "neighbor 1 does not have this as neighbor 0" 
			    << endl;}
	CGAL_triangulation_assertion(false); return false;
      }

      break;
      }

    case 2:
      {
      if ( vertex(0) == NULL || vertex(1) == NULL || vertex(2) == NULL ) {
	if (verbose) { cerr << "vertex 0, 1, or 2 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( vertex(3) != NULL ) {
	if (verbose) { cerr << "vertex 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( neighbor(3) != NULL ) {
	if (verbose) { cerr << "neighbor 3 != NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }

      int in;
      Cell* n;
      for(int i = 0; i < 3; i++) {
	n = neighbor(i);
	if ( n == NULL ) {
	  if (verbose) { cerr << "neighbor " << i << " NULL" << endl;}
	  CGAL_triangulation_assertion(false); return false;
	}
	if ( ! n->has_vertex(vertex(cw(i)),in ) ) {
	  if (verbose) { cerr << "vertex " << cw(i) 
			      << " not vertex of neighbor " << i << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	in = cw(in); 
	if ( this != n->neighbor(in) ) {
	  if (verbose) { cerr << "neighbor " << i
			      << " does not have this as neighbor " << in << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
	if ( vertex(ccw(i)) != n->vertex(cw(in)) ) {
	  if (verbose) { cerr << "vertex " << ccw(i)
			      << " is not vertex " << cw(in) 
			      << " of neighbor " << i << endl;}
	  CGAL_triangulation_assertion(false); return false;
	}
      }
      break;
      }

    case 3:
      {
//     for(int i = 0; i < 4; i++) {
//       if ( vertex(i) == NULL ) return false;
//     }

//     for(int i = 0; i < 4; i++) {
//       Cell* n = neighbor(i);
//       if ( n == NULL ) return false;
            
//       int in;
//       if ( ! n->has_neighbor(this,in) ) return false;

//       int j1n,j2n,j3n;
//       if ( ! n->has_vertex(vertex((i+1)%4),j1n) ) return false;
//       if ( ! n->has_vertex(vertex((i+2)%4),j2n) ) return false;
//       if ( ! n->has_vertex(vertex((i+3)%4),j3n) ) return false;

//       // tests whether the orientations of this and n are consistent
//       switch ( j1n ) {
//       case (in+1)%4 :
// 	if ( ( ! j2n == (in+3)%4 ) || ( ! j3n == (in+2)%4 ) ) return false;
//       case (in+2)%4 :
// 	if ( ( ! j2n == (in+1)%4 ) || ( ! j3n == (in+3)%4 ) ) return false;
//       case (in+3)%4 :
// 	if ( ( ! j2n == (in+2)%4 ) || ( ! j3n == (in+1)%4 ) ) return false;
//       };
//     }
      }
    }
    return true;
  }

private:

  // to maintain the list of cells
  Cell* _previous_cell;
  Cell* _next_cell;
  
  inline
  void add_list(Tds& tds)
  {
    this->_next_cell = tds.list_of_cells()._next_cell;
    tds.list_of_cells()._next_cell=this;
    this->_next_cell->_previous_cell = this;
    this->_previous_cell = tds.past_end_cell();
  }

  // only for dimension 2

  inline int ccw(int i) const
  {
    return (i+1) % 3;
  }
    
  inline int cw(int i) const
  {
    return (i+2) % 3;
  }
 
};

#endif CGAL_TRIANGULATION_DS_CELL_H
