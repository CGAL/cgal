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
// file          : include/CGAL/Triangulation_ds_cell_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// cell of a combinatorial triangulation of any dimension <=3
// use to store vertices if dimension <=0, edges if dimension 1,
// faces if dimension 2, plain cells if dimension 3

#ifndef CGAL_TRIANGULATION_DS_CELL_3_H
#define CGAL_TRIANGULATION_DS_CELL_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_utils_3.h>

CGAL_BEGIN_NAMESPACE

template < class Vb, class Cb > class Triangulation_ds_vertex_3;

template < class Vb, class Cb > class Triangulation_data_structure_3;

template < class Tds> class Triangulation_ds_cell_iterator_3;
template < class Tds> class Triangulation_ds_facet_iterator_3;
template < class Tds> class Triangulation_ds_edge_iterator_3;
template < class Tds> class Triangulation_ds_vertex_iterator_3;

template < class Vb, class Cb >
class Triangulation_ds_cell_3
  : public Cb
{
    static int ccw(int i)
    {
      return Triangulation_utils_3::ccw(i);
    } 

    static int cw(int i)
    {
      return Triangulation_utils_3::cw(i);
    } 

public:

  typedef Triangulation_data_structure_3<Vb,Cb> Tds;

  friend class Triangulation_data_structure_3<Vb,Cb>;

  friend class Triangulation_ds_cell_iterator_3<Tds>;
  friend class Triangulation_ds_facet_iterator_3<Tds>;
  friend class Triangulation_ds_edge_iterator_3<Tds>;
  friend class Triangulation_ds_vertex_iterator_3<Tds>;

  typedef Triangulation_ds_vertex_3<Vb,Cb> Vertex;
  typedef Triangulation_ds_cell_3<Vb,Cb> Cell;

  Triangulation_ds_cell_3()
    : Cb(), in_conflict_flag(0)
  {}

  Triangulation_ds_cell_3(Cell* c)
    : Cb(*c), in_conflict_flag(0)
  {}
    
  Triangulation_ds_cell_3(Vertex* v0, Vertex* v1, 
			  Vertex* v2, Vertex* v3)
    : Cb(v0,v1,v2,v3), in_conflict_flag(0)
  {}

  Triangulation_ds_cell_3(Vertex* v0, Vertex* v1, 
			  Vertex* v2, Vertex* v3,
			  Cell* n0, Cell* n1, Cell* n2, Cell* n3)
    : Cb(v0,v1,v2,v3,n0,n1,n2,n3), in_conflict_flag(0)
  {}

  // SETTING

  void set_vertex(int i, Vertex* v)
  {
    Cb::set_vertex(i,v);
  }

  void set_neighbor(int i, Cell* n)
  {
    Cb::set_neighbor(i,n);
  }

  void set_vertices(Vertex* v0, Vertex* v1, Vertex* v2, Vertex* v3)
  {
    Cb::set_vertices(v0,v1,v2,v3);
  }

  void set_neighbors(Cell* n0, Cell* n1, Cell* n2, Cell* n3)
  {
    Cb::set_neighbors(n0,n1,n2,n3);
  }

  // VERTEX ACCESS

  Vertex* vertex(int i) const
  {
    return (Vertex*) (Cb::vertex(i));
  } 

  bool has_vertex(const Vertex* v) const
  {
    return Cb::has_vertex(v);
  }
    
  bool has_vertex(const Vertex* v, int & i) const
  {
    return Cb::has_vertex(v,i);
  }
    
  int index(const Vertex* v) const
  {
    return Cb::vertex_index(v);
  }

  // NEIGHBOR ACCESS

  Cell* neighbor(int i) const
  {
    return (Cell*) Cb::neighbor(i);
  }
    
  bool has_neighbor(const Cell* n) const
  {
    return Cb::has_neighbor(n);
  }
    
  bool has_neighbor(const Cell* n, int & i) const
  {
    return Cb::has_neighbor(n,i);
  }
    
  int index(const Cell* n) const
  {
    return Cb::cell_index(n);
  }

  int mirror_index(int i) const
  {
      CGAL_triangulation_precondition ( i>=0 && i<4 );
      return neighbor(i)->index(this);
  }
      
  // CHECKING
    
  Vertex* mirror_vertex(int i) const
  {
      return neighbor(i)->vertex(mirror_index(i));
  }

  bool is_valid(int dim = 3, bool verbose = false, int level = 0) const;

private:

  // to maintain the list of cells
  Cell* _previous_cell;
  Cell* _next_cell;
  int in_conflict_flag;
  
  void set_in_conflict_flag(int f)
  {
    in_conflict_flag = f;
  }

  int get_in_conflict_flag() const
  {
    return in_conflict_flag;
  }

  void error_orient( Cell * , int i ) const
  {
    std::cerr << " pb orientation with neighbor " << i << std::endl;
  }

  void error_neighbor( Cell* , int , int ) const
  {
    std::cerr << "neighbor of c has not c as neighbor" << std::endl;
  }
};

template <class Vb, class Cb >
bool
Triangulation_ds_cell_3<Vb,Cb>::is_valid
(int dim, bool verbose, int level) const
{
    if ( ! Cb::is_valid(verbose,level) )
	return false;

    switch (dim) {
    case -2:
    case -1:
    {
      if ( vertex(0) == NULL ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      if ( vertex(1) != NULL || vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose)
	    std::cerr << "vertex 1,2 or 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(0) != NULL || neighbor(1) != NULL ||
	   neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose)
	    std::cerr << "one neighbor != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      break;
    }

    case 0:
      {
      if ( vertex(0) == NULL ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      if ( neighbor (0) == NULL ) {
	if (verbose)
	    std::cerr << "neighbor 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( vertex(1) != NULL || vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose)
	    std::cerr << "vertex 1, 2 or 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(1) != NULL ||
	   neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose)
	    std::cerr << "neighbor 1, 2 or 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      if ( ! neighbor(0)->has_vertex(vertex(0)) ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have vertex 0" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
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
	if (verbose)
	    std::cerr << "vertex 0 or 1 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
      if ( n0 == NULL || n1 == NULL ) {
	if (verbose)
	    std::cerr << "neighbor 0 or 1 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( vertex(2) != NULL || vertex(3) != NULL ) {
	if (verbose)
	    std::cerr << "vertex 2 or 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(2) != NULL || neighbor(3) != NULL ) {
	if (verbose)
	    std::cerr << "neighbor 2 or 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      if ( v0 !=  n1->vertex(1) ) {
	if (verbose)
	    std::cerr << "neighbor 1 does not have vertex 0 as vertex 1"
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( v1 != n0->vertex(0) ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have vertex 1 as vertex 0" 
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      
      if ( this != n0->neighbor(1) ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have this as neighbor 1" 
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( this != n1->neighbor(0) ) {
	if (verbose)
	    std::cerr << "neighbor 1 does not have this as neighbor 0" 
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      break;
      }

    case 2:
      {
      if ( vertex(0) == NULL || vertex(1) == NULL || vertex(2) == NULL ) {
	if (verbose)
	    std::cerr << "vertex 0, 1, or 2 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
      vertex(2)->is_valid(verbose,level);
      if ( vertex(3) != NULL ) {
	if (verbose)
	    std::cerr << "vertex 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(3) != NULL ) {
	if (verbose)
	    std::cerr << "neighbor 3 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }

      int in;
      Cell* n;
      for(int i = 0; i < 3; i++) {
	n = neighbor(i);
	if ( n == NULL ) {
	  if (verbose)
	      std::cerr << "neighbor " << i << " NULL" << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	if ( ! n->has_vertex(vertex(cw(i)),in ) ) {
	  if (verbose)
	      std::cerr << "vertex " << cw(i) 
		        << " not vertex of neighbor " << i << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	in = cw(in); 
	if ( this != n->neighbor(in) ) {
	  if (verbose)
	      std::cerr << "neighbor " << i
		        << " does not have this as neighbor " 
		        << in << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
	if ( vertex(ccw(i)) != n->vertex(cw(in)) ) {
	  if (verbose)
	      std::cerr << "vertex " << ccw(i)
		        << " is not vertex " << cw(in) 
		        << " of neighbor " << i << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
      }
      break;
      }

    case 3:
      {
	int i;
	for(i = 0; i < 4; i++) {
	  if ( vertex(i) == NULL ) {
	    if (verbose)
		std::cerr << "vertex " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  vertex(i)->is_valid(verbose,level);
	}

	for(i = 0; i < 4; i++) {
	  Cell* n = neighbor(i);
	  if ( n == NULL ) {
	    if (verbose)
	      std::cerr << "neighbor " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }

	  int in;
	  if ( ! n->has_neighbor(this,in) ) {
	    if (verbose)
	      error_neighbor(n,i,in); 
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  
	  int j1n,j2n,j3n;
	  if ( ! n->has_vertex(vertex((i+1)&3),j1n) ) {
	    if (verbose) { std::cerr << "vertex " << ((i+1)&3)
				     << " not vertex of neighbor " 
				     << i << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  if ( ! n->has_vertex(vertex((i+2)&3),j2n) ) {
	    if (verbose) { std::cerr << "vertex " << ((i+2)&3)
				     << " not vertex of neighbor " 
				     << i << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  if ( ! n->has_vertex(vertex((i+3)&3),j3n) ) {
	    if (verbose) { std::cerr << "vertex " << ((i+3)&3)
				     << " not vertex of neighbor "
				     << i << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  
	  if ( in+j1n+j2n+j3n != 6) {
	    if (verbose) { std::cerr << "sum of the indices != 6 " 
				     << std::endl; }
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  
	  // tests whether the orientations of this and n are consistent
	  if ( ((i+in)&1) == 0 ) { // i and in have the same parity
	    if ( j1n == ((in+1)&3) ) {
	      if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+2)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	    if ( j1n == ((in+2)&3) ) {
	      if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+3)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	    if ( j1n == ((in+3)&3) ) {
	      if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+1)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	  }
	  else { // i and in do not have the same parity
	    if ( j1n == ((in+1)&3) ) {
	      if ( ( j2n != ((in+2)&3) ) || ( j3n != ((in+3)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	    if ( j1n == ((in+2)&3) ) {
	      if ( ( j2n != ((in+3)&3) ) || ( j3n != ((in+1)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	    if ( j1n == ((in+3)&3) ) {
	      if ( ( j2n != ((in+1)&3) ) || ( j3n != ((in+2)&3) ) ) {
		if (verbose) { 
		  error_orient(n,i);
		}
		CGAL_triangulation_assertion(false);
		return false;
	      }
	    }
	  }
	} // end looking at neighbors
      }// end case dim 3
    } // end switch
    return true;
} // end is_valid

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DS_CELL_3_H
