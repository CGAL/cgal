// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_COUNT_ALPHA_H
#define CGAL_COUNT_ALPHA_H

#include <list>


template <class AS>
void
count_faces(const AS &A, bool verbose)
{
  typedef AS         Alpha_shape_3;


  typedef typename AS::Alpha_shape_cells_iterator
                                      Alpha_shape_cells_iterator;
  typedef typename AS::Alpha_shape_vertices_iterator
                                      Alpha_shape_vertices_iterator;
  typedef typename AS::Alpha_shape_facets_iterator
                                      Alpha_shape_facets_iterator;
  typedef typename AS::size_type   size_type;


  typedef typename AS::Cell_handle   Cell_handle;
  typedef typename AS::Facet         Facet;
  typedef typename AS::Edge          Edge;
  typedef typename AS::Vertex_handle Vertex_handle;
  typedef std::list<Cell_handle>     Cell_list;
  typedef std::list<Facet>           Facet_list;
  typedef std::list<Edge>            Edge_list;
  typedef std::list<Vertex_handle>   Vertex_list;

  Alpha_shape_cells_iterator cit=A.alpha_shape_cells_begin();
  size_type count_cells=0;
  for ( ; cit != A.alpha_shape_cells_end() ; ++cit) ++count_cells;

  Cell_list exterior_cells;
  Cell_list interior_cells;
  A.get_alpha_shape_cells( std::back_inserter(exterior_cells),
                           Alpha_shape_3::EXTERIOR);
  A.get_alpha_shape_cells( std:: back_inserter(interior_cells),
                           Alpha_shape_3::INTERIOR);
  assert(count_cells == interior_cells.size());
  assert(interior_cells.size() + exterior_cells.size() ==
         A.number_of_finite_cells());


  size_type count_facets = 0;
  Facet_list exterior_facets;
  Facet_list singular_facets;
  Facet_list regular_facets;
  Facet_list interior_facets;
  Alpha_shape_facets_iterator face_iterator = A.alpha_shape_facets_begin();
  for (;face_iterator!=A.alpha_shape_facets_end();face_iterator++)
    count_facets++;
  A.get_alpha_shape_facets(std::back_inserter(exterior_facets),
                           Alpha_shape_3::EXTERIOR);
  A.get_alpha_shape_facets(std::back_inserter(singular_facets),
                           Alpha_shape_3::SINGULAR);
  A.get_alpha_shape_facets(std::back_inserter(regular_facets),
                           Alpha_shape_3::REGULAR);
  A.get_alpha_shape_facets(std::back_inserter(interior_facets),
                           Alpha_shape_3::INTERIOR);
  size_type count_exterior_facets = exterior_facets.size();
  size_type count_singular_facets = singular_facets.size();
  size_type count_regular_facets  = regular_facets.size();
  size_type count_interior_facets = interior_facets.size();

  Edge_list exterior_edges;
  Edge_list singular_edges;
  Edge_list regular_edges;
  Edge_list interior_edges;
  A.get_alpha_shape_edges(std::back_inserter(exterior_edges),
                          Alpha_shape_3::EXTERIOR);
  A.get_alpha_shape_edges(std::back_inserter(singular_edges),
                          Alpha_shape_3::SINGULAR);
  A.get_alpha_shape_edges(std::back_inserter(regular_edges),
                          Alpha_shape_3::REGULAR);
  A.get_alpha_shape_edges(std::back_inserter(interior_edges),
                          Alpha_shape_3::INTERIOR);
  size_type count_exterior_edges = exterior_edges.size();
  size_type count_singular_edges = singular_edges.size();
  size_type count_regular_edges  = regular_edges.size();
  size_type count_interior_edges = interior_edges.size();


  size_type count_vertices = 0;
  Alpha_shape_vertices_iterator vit=A.alpha_shape_vertices_begin();
  for ( ; vit != A.alpha_shape_vertices_end() ; ++vit) ++count_vertices;

  Vertex_list exterior_vertices;
  Vertex_list singular_vertices;
  Vertex_list regular_vertices;
  Vertex_list interior_vertices;
  A.get_alpha_shape_vertices(std::back_inserter(exterior_vertices),
                             Alpha_shape_3::EXTERIOR);
  A.get_alpha_shape_vertices(std::back_inserter(singular_vertices),
                             Alpha_shape_3::SINGULAR);
  A.get_alpha_shape_vertices(std::back_inserter(regular_vertices),
                             Alpha_shape_3::REGULAR);
  A.get_alpha_shape_vertices(std::back_inserter(interior_vertices),
                             Alpha_shape_3::INTERIOR);
  size_type count_exterior_vertices = exterior_vertices.size();
  size_type count_singular_vertices = singular_vertices.size();
  size_type count_regular_vertices  = regular_vertices.size();
  size_type count_interior_vertices = interior_vertices.size();


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
              << count_exterior_vertices  << "\t"
              << count_singular_vertices << "\t"
              << count_regular_vertices << "\t"
              << count_interior_vertices << "\t"
              << count_vertices << std::endl;
    std::cerr<< std::endl;
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

   assert(count_singular_vertices == 0);
   assert(count_vertices == count_regular_vertices );
   assert(  count_interior_vertices
           + count_regular_vertices
           + count_singular_vertices
           + count_exterior_vertices == A.number_of_vertices());
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

    assert( count_vertices ==  count_regular_vertices
                             + count_singular_vertices);
    assert(  count_exterior_vertices
           + count_regular_vertices
           + count_singular_vertices
           + count_interior_vertices == A.number_of_vertices());
  }

  if (count_cells >= 1 && A.get_mode()== Alpha_shape_3::REGULARIZED) {
     //this relation might not be valid for any alpha_shape
     // if connected components are touching
     // through an edge or a vertex
     assert(count_regular_facets == 2*count_regular_vertices - 4*A.number_of_solid_components());
   }

}

// forward declaration
template <class AS>
void show_alpha_status(AS&, const typename AS::Alpha_status&);

template <class AS>
void
test_filtration(AS &A, bool verbose)
{
  typename std::list<CGAL::Object> filtration;

  A.filtration(std::back_inserter(filtration));
  typename AS::size_type count_vertices = 0;
  typename AS::size_type count_edges = 0;
  typename AS::size_type count_facets = 0;
  typename AS::size_type count_cells = 0;



    typename std::list<CGAL::Object>::iterator filtre_it = filtration.begin();
    typename AS::Cell_handle cell;
    typename AS::Facet      facet;
    typename AS::Edge       edge;
    typename AS::Vertex_handle  vertex;
    typename AS::Alpha_status  as;
    typename AS::NT alpha;
    if(verbose) {
      std::cerr << std::endl;
      std::cerr << "Analyze filtration " << std::endl;
    }
    for (; filtre_it != filtration.end(); filtre_it++) {
      if(assign(vertex, *filtre_it)) {
        as = *(vertex->get_alpha_status());
        if(verbose) std::cerr << "Vertex" << "\t";
        if(verbose)show_alpha_status(A,as);
        count_vertices++;
      }
      if(assign(edge, *filtre_it)) {
        // could be done with Edge_alpha_map in GENERAL mode
        A.compute_edge_status(edge.first,edge.second,edge.third, as);
        if(verbose) std::cerr << "Edge" << "\t";
        if(verbose) show_alpha_status(A,as);
        count_edges++;
      }
      if(assign(facet, *filtre_it)) {
        as = *(facet.first->get_facet_status(facet.second));
        if(verbose) std::cerr << "Facet" << "\t";
        if(verbose) show_alpha_status(A,as);
        count_facets++;
      }
      if(assign(cell, *filtre_it)) {
        alpha  = cell->get_alpha();
        if(verbose) std::cerr << "Cell" << "\t";
        if(verbose) std::cerr << alpha << std::endl;
        count_cells++;
      }
    }
    if(verbose) {
      std::cerr << "vertices \t" << count_vertices << "\t"
                << A.number_of_vertices() << std::endl;
      std::cerr << "edges \t" << count_edges << "\t"
                << A.number_of_finite_edges() << std::endl;
      std::cerr << "facets \t" << count_facets << "\t"
                << A.number_of_finite_facets() << std::endl;
      std::cerr << "cellss \t" << count_cells << "\t"
                << A.number_of_finite_cells() << std::endl;
    }
    assert(count_vertices == A.number_of_vertices());
    assert(count_edges == A.number_of_finite_edges());
    assert(count_facets == A.number_of_finite_facets());
    assert(count_cells ==  A.number_of_finite_cells());

    filtration.clear();
    std::list<typename AS::FT> alpha_values;

    A.filtration_with_alpha_values(
      CGAL::dispatch_output<CGAL::Object,typename AS::FT>(
        std::back_inserter(filtration),
        std::back_inserter(alpha_values)
      )
    );

    assert( filtration.size() == alpha_values.size() );
    typename std::list<typename AS::FT>::iterator av_it=alpha_values.begin();
    typename AS::FT previous_alpha_value=*av_it++;
    for(;av_it!=alpha_values.end();++av_it)
    {
      assert( previous_alpha_value <= *av_it );
      previous_alpha_value=*av_it;
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
void show_alpha_values(AS& A)
{
  typedef typename AS::Alpha_iterator Alpha_iterator;

  std::cerr <<"Alpha spectrum \n";
  for (Alpha_iterator alpha_it = A.alpha_begin();
       alpha_it!=A.alpha_end(); alpha_it++){
    std::cerr<< *alpha_it <<std::endl;
  }

}

template <class AS>
void show_alpha_status(AS& A, const typename AS::Alpha_status&  as)
{
  if (A.get_mode() == AS::REGULARIZED || !(as.is_Gabriel())){
    std::cerr << "*****" << "\t" ;
  } else{
    std::cerr << as.alpha_min() << "\t";
  }
  std::cerr << as.alpha_mid() << "\t";
  if (as.is_on_chull())  std::cerr << "*****" << std::endl;
  else  std::cerr << as.alpha_max() << std::endl;
  return;
}


#endif //  CGAL_COUNT_ALPHA_H
