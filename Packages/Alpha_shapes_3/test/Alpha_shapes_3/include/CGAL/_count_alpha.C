// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
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
// $Revision$ 
// $Date$
// $Name$
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_COUNT_ALPHA_C
#define CGAL_COUNT_ALPHA_C


template <class AS>
void 
count_faces(const AS &A, bool verbose)
{
  typedef AS         Alpha_shape_3;

  typedef typename AS::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename AS::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename AS::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename AS::Finite_vertices_iterator Finite_vertices_iterator;
  
  typedef typename AS::Alpha_shape_cells_iterator 
                                      Alpha_shape_cells_iterator;
  typedef typename AS::Alpha_shape_vertices_iterator 
                                      Alpha_shape_vertices_iterator;
  typedef typename AS::Alpha_shape_facets_iterator 
                                      Alpha_shape_facets_iterator; 
  typedef typename AS::size_type   size_type;

  Alpha_shape_cells_iterator cit=A.alpha_shape_cells_begin();
  size_type count_cells=0;
  for ( ; cit != A.alpha_shape_cells_end() ; ++cit) ++count_cells;
  Finite_cells_iterator fcit=A.finite_cells_begin();
  size_type count_cint=0;
  for (; fcit != A.finite_cells_end(); ++fcit) 
    if ( A.classify(fcit) ==  Alpha_shape_3::INTERIOR) ++count_cint;
  assert(count_cells == count_cint);

  size_type count_exterior_facets = 0;
  size_type count_singular_facets = 0;
  size_type count_regular_facets  = 0;
  size_type count_interior_facets = 0;
  size_type count_facets = 0;
  Alpha_shape_facets_iterator face_iterator = A.alpha_shape_facets_begin();
  for (;face_iterator!=A.alpha_shape_facets_end();face_iterator++) 
    count_facets++;
  Finite_facets_iterator fit = A.finite_facets_begin();
  for ( ;  fit != A.finite_facets_end(); ++fit) {
    switch(A.classify(*fit) ) {
    case Alpha_shape_3::EXTERIOR : ++count_exterior_facets; break;
    case Alpha_shape_3::SINGULAR : ++count_singular_facets;  break;
    case Alpha_shape_3::REGULAR  : ++count_regular_facets;  break;
    case Alpha_shape_3::SUPER_REGULAR : break;
    case Alpha_shape_3::INTERIOR : ++count_interior_facets;  break;
    }
  }
  
  size_type count_exterior_edges = 0;
  size_type count_singular_edges = 0;
  size_type count_regular_edges = 0;
  size_type count_interior_edges = 0;
  for (Finite_edges_iterator e_iterator = A.finite_edges_begin();
       e_iterator!=A.finite_edges_end();e_iterator++){
       switch(A.classify(*e_iterator)) {
       case Alpha_shape_3::EXTERIOR : ++count_exterior_edges; break;
       case Alpha_shape_3::SINGULAR : ++count_singular_edges;  break;
       case Alpha_shape_3::REGULAR  : ++count_regular_edges;  break;
       case Alpha_shape_3::SUPER_REGULAR : break;
       case Alpha_shape_3::INTERIOR : ++count_interior_edges;  break;
    }
  }

  size_type count_regular_vertices = 0;
  size_type count_super_regular_vertices = 0;
  size_type count_singular_vertices = 0;
  size_type count_exterior_vertices = 0;
  size_type count_interior_vertices = 0;
  size_type count_vertices = 0;
  Alpha_shape_vertices_iterator vit=A.alpha_shape_vertices_begin();
  for ( ; vit != A.alpha_shape_vertices_end() ; ++vit) ++count_vertices;
  Finite_vertices_iterator fvit = A.finite_vertices_begin();
  for ( ; fvit != A.finite_vertices_end(); ++fvit) {
    switch(A.classify(fvit)) {
    case Alpha_shape_3::EXTERIOR : ++count_exterior_vertices; break;
    case Alpha_shape_3::SINGULAR : ++count_singular_vertices;  break;
    case Alpha_shape_3::REGULAR  : ++count_regular_vertices;  break;
    case Alpha_shape_3::SUPER_REGULAR  : 
                              ++count_super_regular_vertices;  break;
    case Alpha_shape_3::INTERIOR : ++count_interior_vertices;  break;
    }
  }

  if (verbose) {
    std::cerr << "facets " << "\t" << "\t"
	      << count_exterior_facets << "\t"
	      << count_singular_facets << "\t" 
	      << count_regular_facets << "\t"
	      << count_interior_facets << "\t"
	      << count_facets << std::endl;

    std::cerr << "edges " << "\t" << "\t"
	      << count_exterior_edges << "\t"
	      << count_singular_edges << "\t"
	      << count_regular_edges << "\t"
	      << count_interior_edges << std::endl;

    std::cerr << "vertices "<< "\t"
	      << count_singular_vertices  << "\t"
	      << count_regular_vertices << "\t" 
	      << count_super_regular_vertices << "\t" 
	      << count_interior_vertices << "\t"
	      << count_vertices << std::endl;
  }


  if (A.get_mode() == Alpha_shape_3::REGULARIZED){
    assert( count_facets == count_regular_facets );
    assert (count_singular_facets == 0);
    assert (count_exterior_facets 
	    + count_regular_facets 
	    + count_interior_facets == A.number_of_finite_facets());

   assert(count_singular_edges == 0); 
   assert(  count_interior_edges
	   + count_regular_edges
	   + count_exterior_edges == A.number_of_finite_edges());

  assert(count_regular_vertices == 0);
  assert(count_exterior_vertices == 0);
  assert(count_vertices == count_super_regular_vertices );
  assert(  count_interior_vertices 
	   + count_super_regular_vertices 
	   + count_singular_vertices  == A.number_of_vertices());
  }

  if( A.get_mode() == Alpha_shape_3::GENERAL) {
    assert( count_facets == count_regular_facets + count_singular_facets);
    assert (  count_exterior_facets 
	    + count_singular_facets
	    + count_regular_facets 
	    + count_interior_facets == A.number_of_finite_facets());

    assert(    count_interior_edges
	     + count_regular_edges
	     + count_singular_edges 
	     + count_exterior_edges == A.number_of_finite_edges());

    assert(count_exterior_vertices == 0);
    assert( count_vertices == 
	    count_super_regular_vertices
	    + count_regular_vertices
	    + count_singular_vertices);
    assert(  count_super_regular_vertices
	   + count_regular_vertices
	   + count_singular_vertices
	   + count_interior_vertices == A.number_of_vertices());
  }

  int ncc = A.number_of_solid_components();

  if (count_cells >= 1) {
   //this relation might not be valid for any alpha_shape
   // if connected components are touching
   // through an edge or a vertex
   assert(count_regular_facets == 2*count_super_regular_vertices - 4*ncc);
 }
}


//TO DEBUG
template <class AS>
void
show_triangulation(AS& A)
{
  std::cerr << "number of finite cells " << A.number_of_finite_cells() 
	    << std::endl;
  std::cerr << "number of finite faces " << A.number_of_finite_facets() 
	    << std::endl;
  std::cerr << "number of finite edges " << A.number_of_finite_edges() 
	    << std::endl;

  typedef typename AS::Finite_cells_iterator Finite_cells_iterator;
  Finite_cells_iterator cit = A.finite_cells_begin();
  for (  Finite_cells_iterator cit = A.finite_cells_begin();
	 cit != A.finite_cells_end(); ++cit) {
    std::cerr << "cell " <<  "alpha " << cit->get_alpha() << std::endl;
    std::cerr << cit->vertex(0)->point() << std::endl
	      << cit->vertex(1)->point() << std::endl
	      << cit->vertex(2)->point() << std::endl
	      << cit->vertex(3)->point() << std::endl;
  }
}

template <class AS>
void show_alpha_values(AS& A )
{
  typedef typename AS::Alpha_iterator Alpha_iterator;

  std::cerr <<"Alpha spectrum \n";
  for (Alpha_iterator alpha_it = A.alpha_begin();
       alpha_it!=A.alpha_end(); alpha_it++){
    std::cerr<< *alpha_it <<std::endl;
  }

}


#endif //  CGAL_COUNT_ALPHA_C
