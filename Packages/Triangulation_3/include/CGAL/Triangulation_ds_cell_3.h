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
// file          : include/CGAL/Triangulation_ds_cell_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================
//
// cell of a combinatorial triangulation of any dimension <=3
// use to store vertices if dimension <=0, edges if dimension 1,
// faces if dimension 2, plain cells if dimension 3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_CELL_3_H
#define CGAL_TRIANGULATION_DS_CELL_3_H

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_utils_3.h>

template < class Vb, class Cb >
class CGAL_Triangulation_ds_vertex_3;

template < class Vb, class Cb >
class CGAL_Triangulation_data_structure_3;

//template < class Tds>
//class CGAL_Triangulation_ds_iterator_base_3;
template < class Tds>
class CGAL_Triangulation_ds_cell_iterator_3;
template < class Tds>
class CGAL_Triangulation_ds_facet_iterator_3;
template < class Tds>
class CGAL_Triangulation_ds_edge_iterator_3;
template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator_3;

template < class Vb, class Cb >
class CGAL_Triangulation_ds_cell_3
  : public Cb,
    public CGAL_Triangulation_utils_3
{

  friend class CGAL_Triangulation_data_structure_3<Vb,Cb>;

  //  friend class CGAL_Triangulation_ds_iterator_base_3
  //  <CGAL_Triangulation_data_structure_3<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_cell_iterator_3
  <CGAL_Triangulation_data_structure_3<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_facet_iterator_3
  <CGAL_Triangulation_data_structure_3<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_edge_iterator_3
  <CGAL_Triangulation_data_structure_3<Vb,Cb> >;
  friend class CGAL_Triangulation_ds_vertex_iterator_3
  <CGAL_Triangulation_data_structure_3<Vb,Cb> >;
  

public:

  //  typedef typename Vb::Point Point;

  typedef CGAL_Triangulation_data_structure_3<Vb,Cb> Tds;
  typedef CGAL_Triangulation_ds_vertex_3<Vb,Cb> Vertex;
  //  typedef typename CGAL_Triangulation_data_structure_3<Vb,Cb>::Facet Facet;
  typedef CGAL_Triangulation_ds_cell_3<Vb,Cb> Cell;

  // CONSTRUCTORS

  // used only for initializing _list_of_cells by the constructors of
  // CGAL_Triangulation_data_structure_3
  // private ?
  inline
  CGAL_Triangulation_ds_cell_3()
    : Cb(), _previous_cell(this), _next_cell(this)
  {}

  inline
  CGAL_Triangulation_ds_cell_3(Tds & tds)
    : Cb()
    // builds a new cell of Triangulation_data_structure_3 and maintains the list of cells
  { add_list(tds); }

  CGAL_Triangulation_ds_cell_3(Tds & tds, Cell* c)
    : Cb(c->vertex(0),c->vertex(1),c->vertex(2),c->vertex(3),
	 c->neighbor(0),c->neighbor(1),c->neighbor(2),c->neighbor(3))
  { add_list(tds); }
    
  CGAL_Triangulation_ds_cell_3(Tds & tds,
			     Vertex* v0, Vertex* v1, 
			     Vertex* v2, Vertex* v3)
    :  Cb(v0,v1,v2,v3)
  { add_list(tds); }
    
  CGAL_Triangulation_ds_cell_3(Tds & tds,
			     Vertex* v0, Vertex* v1, 
			     Vertex* v2, Vertex* v3,
			     Cell* n0, Cell* n1, Cell* n2, Cell* n3)
    :  Cb(v0,v1,v2,v3,n0,n1,n2,n3)
  { add_list(tds); }


  // DESTRUCTOR

  ~CGAL_Triangulation_ds_cell_3()
    {
      _previous_cell->_next_cell = _next_cell;
      _next_cell->_previous_cell = _previous_cell;
      // automatically calls the destructor of the cell base ?
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
  bool has_vertex(const Vertex* v, int & i) const
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
  bool has_neighbor(const Cell* n, int & i) const
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
    if ( ! Cb::is_valid(verbose, true) ) return false;

    switch (dim) {
      
    case -2:

    case -1:
    {
      if ( vertex(0) == NULL ) {
	if (verbose) { cerr << "vertex 0 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      vertex(0)->is_valid(verbose,level);
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
      vertex(0)->is_valid(verbose,level);
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
      Vertex* v0 = vertex(0); 
      Vertex* v1 = vertex(1);
      Cell* n0 = neighbor(0); 
      Cell* n1 = neighbor(1);

      if ( v0 == NULL || v1 == NULL ) {
	if (verbose) { cerr << "vertex 0 or 1 NULL" << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
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

      if ( v0 !=  n1->vertex(1) ) {
	if (verbose) { cerr << 
		       "neighbor 1 does not have vertex 0 as vertex 1" 
			    << endl;}
	CGAL_triangulation_assertion(false); return false;
      }
      if ( v1 != n0->vertex(0) ) {
	if (verbose) { cerr << 
			 "neighbor 0 does not have vertex 1 as vertex 0" 
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
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
      vertex(2)->is_valid(verbose,level);
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
	int i;
	for(i = 0; i < 4; i++) {
	  if ( vertex(i) == NULL ) {
	    if (verbose) { 
	      cerr << "vertex " << i << " NULL" << endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  vertex(i)->is_valid(verbose,level);
	}

	for(i = 0; i < 4; i++) {
	  Cell* n = neighbor(i);
	  if ( n == NULL ) {
	    if (verbose) { 
	      cerr << "neighbor " << i << " NULL" << endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }

	  int in;
	  if ( ! n->has_neighbor(this,in) ) {
	    if (verbose) { 
	      error_neighbor(n,i,in); 
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  
	  int j1n,j2n,j3n;
	  if ( ! n->has_vertex(vertex((i+1)&3),j1n) ) {
	    if (verbose) { cerr << "vertex " << (i+1)%4
				<< " not vertex of neighbor " << i << endl; }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  if ( ! n->has_vertex(vertex((i+2)&3),j2n) ) {
	    if (verbose) { cerr << "vertex " << (i+2)%4
				<< " not vertex of neighbor " << i << endl; }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  if ( ! n->has_vertex(vertex((i+3)&3),j3n) ) {
	    if (verbose) { cerr << "vertex " << (i+3)%4
				<< " not vertex of neighbor " << i << endl; }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  
	  if ( in+j1n+j2n+j3n != 6) {
	    if (verbose) { cerr << "sum of the indices != 6 " << endl; }
	    CGAL_triangulation_assertion(false); return false;
	  }
	  
	  // tests whether the orientations of this and n are consistent
	  if ( ((i+in)&1) == 0 ) { // i and in have the same parity
	    if ( j1n == ((in+1)&3) ) {
	      if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+2)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	    if ( j1n == ((in+2)&3) ) {
	      if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+3)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	    if ( j1n == ((in+3)&3) ) {
	      if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+1)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	  }
	  else { // i and in do not have the same parity
	    if ( j1n == ((in+1)&3) ) {
	      if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+3)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	    if ( j1n == ((in+2)&3) ) {
	      if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+1)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	    if ( j1n == ((in+3)&3) ) {
	      if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+2)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	  }
	} // end looking at neighbors
      }// end case dim 3
    } // end switch
    return true;
  } // end is_valid

private:

  // to maintain the list of cells
  Cell* _previous_cell;
  Cell* _next_cell;
  
  inline
  void add_list(Tds & tds)
  {
    this->_next_cell = tds.list_of_cells()._next_cell;
    tds.list_of_cells()._next_cell=this;
    this->_next_cell->_previous_cell = this;
    this->_previous_cell = tds.past_end_cell();
  }

  // only for dimension 2

//   inline int ccw(int i) const
//   {
//     return (i+1) % 3;
//   }
    
//   inline int cw(int i) const
//   {
//     return (i+2) % 3;
//   }
 

  void error_orient( Cell * n, int i) const
  {
//     cerr << this->vertex(0)->point() << ", "
// 	 << this->vertex(1)->point() << ", "
// 	 << this->vertex(2)->point() << ", "
// 	 << this->vertex(3)->point() << endl
// 	 << " pb orientation with neighbor " << i
// 	 << " : " << endl 
// 	 << n->vertex(0)->point() << ", "
// 	 << n->vertex(1)->point() << ", "
// 	 << n->vertex(2)->point() << ", "
// 	 << n->vertex(3)->point() << endl
// 	 << endl;
    cerr << " pb orientation with neighbor " << endl;
  }

  void error_neighbor( Cell* n, int i, int in ) const
  {
//     cerr << "neighbor " << i << endl
// 	 << n->vertex(0)->point() << ", "
// 	 << n->vertex(1)->point() << ", "
// 	 << n->vertex(2)->point() << ", "
// 	 << n->vertex(3)->point() << endl
// 	 << " does not have this " << endl
// 	 << this->vertex(0)->point() << ", "
// 	 << this->vertex(1)->point() << ", "
// 	 << this->vertex(2)->point() << ", "
// 	 << this->vertex(3)->point() << endl
// 	 << " as neighbor " << in << endl;
    cerr << "neighbor of c has not c as neighbor" << endl;
  }

};

#endif CGAL_TRIANGULATION_DS_CELL_3_H
