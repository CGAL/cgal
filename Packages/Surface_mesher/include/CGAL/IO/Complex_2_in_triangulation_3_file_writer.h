// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT


#ifndef _COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
#define _COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H


#include <iomanip>


template < class Tr>
int number_of_facets_on_surface(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit)
    if (fit->first->is_facet_on_surface (fit->second))
      ++result;
  return result;
}


template < class Tr>
void
output_surface_facets_to_off (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // Header.


  os << "OFF \n"
     << T.number_of_vertices() << " " <<
    number_of_facets_on_surface (T) << 
    " " << 0 << "\n";

  os << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }
  
  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	os << "3 ";
	for (int i=0; i<4; i++)
          if (i != (*fit).second)
	    os << V[(*fit).first->vertex(i)] << " ";
	
	os << "\n"; // without color.
      }
}

// only if cells have is_in_domain() method.
template < class Tr>
void
output_oriented_surface_facets_to_off (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // Header.


  os << "OFF \n"
     << T.number_of_vertices() << " " <<
    number_of_facets_on_surface (T) << 
    " " << 0 << "\n";

  os << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os << p.x() << " " << p.y() << " " << p.z() << "\n";
  }
  
  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	typename Tr::Facet f = *fit;
	typename Tr::Facet opposite = T.compute_opposite_facet(f);
	CGAL_assertion (f.first->is_in_domain() !=
			opposite.first->is_in_domain());
	if(!f.first->is_in_domain())
	  f = T.compute_opposite_facet(f);
	os << "3 "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,0))] << " "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,1))] << " "
	   << V[f.first->vertex(T.vertex_triple_index(f.second,2))] << " "	
	   << "\n"; // without color.
      }
}

// only if cells have is_in_domain() method.
template < class Tr>
void
output_surface_facets_to_ghs   (std::ostream& os_points,
				std::ostream& os_faces,
				const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;
  typedef typename Tr::Facet Facet;

  // Header.


  os_points << T.number_of_vertices() << "\n";

  os_points << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    Point p = static_cast<Point>(vit->point());
    os_points << p.x() << " " << p.y() << " " << p.z() << " 0\n";
  }
  
  os_faces << number_of_facets_on_surface(T) << "\n";
  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	Facet f = *fit;
	if(!f.first->is_in_domain())
	  f = T.compute_opposite_facet(f);
	os_faces 
	  << "3 "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,0))] << " "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,1))] << " "
	  << V[f.first->vertex(T.vertex_triple_index(f.second,2))] << " "
	  << "0 0 0 0\n";
      }
}

template < class Tr>
int number_of_facets_in_domain(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit) {
    typename Tr::Cell_handle neighb = fit->first->neighbor (fit->second);
    if ((fit->first->is_in_domain () || neighb->is_in_domain()) &&
	!fit->first->is_facet_on_surface (fit->second))
      ++result;
  }
  return result;
}

template < class Tr>
void
output_interior_facets_to_off (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;

  // Header.


  os << "OFF \n"
     << T.number_of_vertices() << " " <<
    number_of_facets_in_domain (T) << 
    " " << 0 << "\n";

  os << std::setprecision(20);
 
  // Finite vertices coordinates.
  std::map<Vertex_handle, int> V;
  int inum = 0;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) {
    V[vit] = inum++;
    os << vit->point() << "\n";
  }
  
  // Finite facets indices.
  for( Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit){
    typename Tr::Cell_handle neighb = fit->first->neighbor (fit->second);
    if ((fit->first->is_in_domain () || neighb->is_in_domain()) &&
	!fit->first->is_facet_on_surface (fit->second))
      {
	os << "3 ";
	for (int i=0; i<4; i++)
          if (i != (*fit).second)
	    os << V[(*fit).first->vertex(i)] << " ";
	
	os << "\n"; // without color.
      }
  }
}


#endif  // _COMPLEX_2_IN_TRIANGULATION_3_FILE_WRITER_H
