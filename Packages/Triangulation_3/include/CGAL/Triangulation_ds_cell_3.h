// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

// cell of a combinatorial triangulation of any dimension <=3
// use to store vertices if dimension <=0, edges if dimension 1,
// faces if dimension 2, plain cells if dimension 3

#ifndef CGAL_TRIANGULATION_DS_CELL_3_H
#define CGAL_TRIANGULATION_DS_CELL_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_utils_3.h>

CGAL_BEGIN_NAMESPACE

template < class Cb >
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

  typedef typename Cb::Triangulation_data_structure Tds;

public:
  typedef typename Tds::Vertex_handle  Vertex_handle;
  typedef typename Tds::Cell_handle    Cell_handle;
  typedef typename Tds::Vertex         Vertex;
  typedef typename Tds::Cell           Cell;

  Triangulation_ds_cell_3()
    : Cb()
  {
      set_in_conflict_flag(0);
  }

  Triangulation_ds_cell_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3)
  {
      set_in_conflict_flag(0);
  }

  Triangulation_ds_cell_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3,
                          Cell_handle   n0, Cell_handle   n1,
                          Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3)
  {
      set_in_conflict_flag(0);
  }

  int mirror_index(int i) const
  {
      CGAL_triangulation_precondition ( i>=0 && i<4 );
      Cell_handle ni = neighbor(i);
      if (&*ni->neighbor(0) == this) return 0;
      if (&*ni->neighbor(1) == this) return 1;
      if (&*ni->neighbor(2) == this) return 2;
      CGAL_triangulation_assertion(&*ni->neighbor(3) == this);
      return 3;
  }

  Vertex_handle mirror_vertex(int i) const
  {
      return neighbor(i)->vertex(mirror_index(i));
  }

  void set_in_conflict_flag(unsigned char f)
  {
      _in_conflict_flag = f;
  }

  unsigned char get_in_conflict_flag() const
  {
      return _in_conflict_flag;
  }

  bool is_valid(int dim = 3, bool verbose = false, int level = 0) const;

private:

  unsigned char _in_conflict_flag;

  void error_orient(Cell_handle, int i) const
  {
    std::cerr << " pb orientation with neighbor " << i << std::endl;
  }
};

template < class Cb >
bool
Triangulation_ds_cell_3<Cb>::is_valid(int dim, bool verbose, int level) const
{
    if ( ! Cb::is_valid(verbose,level) )
	return false;

    switch (dim) {
    case -2:
    case -1:
    {
      if ( vertex(0) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      if ( vertex(1) != Vertex_handle() || vertex(2) != Vertex_handle()) {
	if (verbose)
	    std::cerr << "vertex 1 or 2 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(0) != Cell_handle() ||
	   neighbor(1) != Cell_handle() ||
	   neighbor(2) != Cell_handle()) {
	if (verbose)
	    std::cerr << "one neighbor != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      break;
    }

    case 0:
      {
      if ( vertex(0) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      if ( neighbor (0) == Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 0 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( vertex(1) != Vertex_handle() || vertex(2) != Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 1 or 2 != NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( neighbor(1) != Cell_handle() || neighbor(2) != Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 1 or 2 != NULL" << std::endl;
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
      Vertex_handle v0 = vertex(0); 
      Vertex_handle v1 = vertex(1);
      Cell_handle n0 = neighbor(0); 
      Cell_handle n1 = neighbor(1);

      if ( v0 == Vertex_handle() || v1 == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0 or 1 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
      if ( n0 == Cell_handle() || n1 == Cell_handle() ) {
	if (verbose)
	    std::cerr << "neighbor 0 or 1 NULL" << std::endl;
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
      
      if ( &*n0->neighbor(1) != this ) {
	if (verbose)
	    std::cerr << "neighbor 0 does not have this as neighbor 1" 
		      << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      if ( &*n1->neighbor(0) != this ) {
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
      if ( vertex(0) == Vertex_handle() ||
	   vertex(1) == Vertex_handle() ||
	   vertex(2) == Vertex_handle() ) {
	if (verbose)
	    std::cerr << "vertex 0, 1, or 2 NULL" << std::endl;
	CGAL_triangulation_assertion(false);
	return false;
      }
      vertex(0)->is_valid(verbose,level);
      vertex(1)->is_valid(verbose,level);
      vertex(2)->is_valid(verbose,level);
      int in;
      Cell_handle n;
      for(int i = 0; i < 3; i++) {
	n = neighbor(i);
	if ( n == Cell_handle() ) {
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
	if ( &*n->neighbor(in) != this ) {
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
	  if ( vertex(i) == Vertex_handle() ) {
	    if (verbose)
		std::cerr << "vertex " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	  vertex(i)->is_valid(verbose,level);
	}

	for(i = 0; i < 4; i++) {
	  Cell_handle n = neighbor(i);
	  if ( n == Cell_handle() ) {
	    if (verbose)
	      std::cerr << "neighbor " << i << " NULL" << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }

	  int in = 5;
	  // if ( ! n->has_neighbor(handle(), in) ) {
          if ( &*n->neighbor(0) == this) in = 0;
          if ( &*n->neighbor(1) == this) in = 1;
          if ( &*n->neighbor(2) == this) in = 2;
          if ( &*n->neighbor(3) == this) in = 3;
          if (in == 5) {
	    if (verbose)
              std::cerr << "neighbor of c has not c as neighbor" << std::endl;
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
